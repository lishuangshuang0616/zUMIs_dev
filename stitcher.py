#!/usr/bin/env python
# V 2.0
# anton jm larsson anton.larsson@ki.se
import argparse
import pysam
import warnings
import numpy as np
import portion as P
import itertools
import time
import os
import json
from scipy.special import logsumexp
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
__version__ = '2.0'
nucleotides = ['A', 'T', 'C', 'G']
nuc_dict = {'A':0, 'T':1, 'C':2, 'G':3, 'N': 4}
np.seterr(divide='ignore')
ll_this_correct = {i:np.log(1-10**(-float(i)/10)) for i in range(1,94)}
ln_3 = np.log(3)
ll_other_correct = {i:-(float(i)*np.log(10))/10 - ln_3 for i in range(1,94)}
ll_N = -np.log(4)
def make_ll_array(e):
    y = np.array([e[0]/3,e[0]/3,e[0]/3,e[0]/3])
    if e[1] != 4:
        y[e[1]] = 1-e[0]
    return np.log10(y)

# taken from https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def intervals_extract(iterable): 
    iterable = sorted(set(iterable)) 
    for _key, group in itertools.groupby(enumerate(iterable), 
    lambda t: t[1] - t[0]): 
        group = list(group) 
        yield [group[0][1], group[-1][1]] 

def interval(t):
    return P.from_data([(True,i[0],i[1], True) for i in t])

def get_time_formatted(time):
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    s = ''.join(['{} day{}, '.format(day, 's'*(1 != day))*(0 != day), 
                 '{} hour{}, '.format(hour,'s'*(1 != hour))*(0 != hour), 
                 '{} minute{}, '.format(minutes,'s'*(1 != minutes))*(0 != minutes), 
                 '{:.2f} second{}, '.format(seconds,'s'*(1 != seconds))*(0 != seconds)])
    s = s[:-2]
    s = s + '.'
    return s

def get_insertions_locs(cigtuples):
    insertion_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 1:
            insertion_locs.extend(range(l, l + c[1]))
            l += c[1]
    return insertion_locs

def get_skipped_tuples(cigtuples, ref_positions):
    skipped_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 3:
            skipped_locs.append((ref_positions[l-1]+1, ref_positions[l]-1))
    return skipped_locs

def using_indexed_assignment(x):
    "https://stackoverflow.com/a/5284703/190597 (Sven Marnach)"
    result = np.empty(len(x), dtype=int)
    temp = x.argsort()
    result[temp] = np.arange(len(x))
    return result

def stitch_reads(read_d, single_end, cell, gene, umi, UMI_tag):
    master_read = {}
    nreads = len(read_d)
    reverse_read1 = []
    exonic_list = [0]*nreads
    intronic_list = [0]*nreads
    seq_list = []
    qual_list = []
    ref_pos_set = set()
    ref_pos_list = []
    for i,read in enumerate(read_d):
        if read.has_tag('GE'):
            exonic = True
        else:
            exonic = False
        if read.has_tag('GI'):
            intronic = True
        else:
            intronic = False
        try:
            Q_list = list(read.query_alignment_qualities)
        except TypeError:
            Q_list = [read.query_alignment_qualities]
        
        seq = read.query_alignment_sequence
        ref_positions_full = read.get_reference_positions(full_length=True)
        if seq is None or ref_positions_full is None:
            continue

        filtered_seq = []
        filtered_quals = []
        ref_positions = []
        for pos, base, q in zip(ref_positions_full, seq, Q_list):
            if pos is None:
                continue
            ref_positions.append(pos)
            filtered_seq.append(base)
            filtered_quals.append(q)

        seq = "".join(filtered_seq)
        Q_list = filtered_quals
        blocks = read.get_blocks()
        skipped_intervals = []
        if blocks is not None and len(blocks) > 1:
            for n in range(len(blocks) - 1):
                s = blocks[n][1]
                e = blocks[n + 1][0] - 1
                if s <= e:
                    skipped_intervals.append((s, e))

        if read.is_read1 and not single_end and read.get_tag(UMI_tag) != '':
            reverse_read1.append(read.is_reverse)
        elif single_end:
            reverse_read1.append(read.is_reverse)

        exonic_list[i] = exonic
        intronic_list[i] = intronic

        seq_list.append(seq)

        qual_list.append(Q_list)

        ref_pos_list.append(ref_positions)

        ref_pos_set.update(ref_positions)


        if len(master_read) == 0:
            master_read['skipped_intervals'] = skipped_intervals
        else:
            master_read['skipped_intervals'].extend(skipped_intervals)

    ref_pos_set_array = np.array(list(ref_pos_set), dtype=np.int64)
    if ref_pos_set_array.size == 0:
        return (False, ':'.join([gene,cell,umi]))

    sorted_ref = np.sort(ref_pos_set_array)

    ll_this_arr = np.full(256, -np.inf, dtype=np.float64)
    ll_other_arr = np.full(256, -np.inf, dtype=np.float64)
    qs = np.arange(0, 94, dtype=np.float64)
    ll_this_arr[:94] = np.log(1 - np.power(10.0, -qs / 10.0))
    ll_other_arr[:94] = -(qs * np.log(10.0)) / 10.0 - np.log(3.0)

    base_to_idx = np.full(256, -1, dtype=np.int8)
    base_to_idx[ord('A')] = 0
    base_to_idx[ord('T')] = 1
    base_to_idx[ord('C')] = 2
    base_to_idx[ord('G')] = 3

    llA = np.zeros(sorted_ref.size, dtype=np.float64)
    llT = np.zeros(sorted_ref.size, dtype=np.float64)
    llC = np.zeros(sorted_ref.size, dtype=np.float64)
    llG = np.zeros(sorted_ref.size, dtype=np.float64)

    for seq, Q_list, ref_positions in zip(seq_list, qual_list, ref_pos_list):
        if not ref_positions:
            continue
        pos_arr = np.asarray(ref_positions, dtype=np.int64)
        cols = np.searchsorted(sorted_ref, pos_arr)

        seq_b = seq.encode('ascii')
        base_arr = np.frombuffer(seq_b, dtype=np.uint8)
        q_arr = np.asarray(Q_list, dtype=np.uint8)

        mask_n = base_arr == ord('N')
        if np.any(mask_n):
            cols_n = cols[mask_n]
            llA[cols_n] += ll_N
            llT[cols_n] += ll_N
            llC[cols_n] += ll_N
            llG[cols_n] += ll_N

        mask_b = ~mask_n
        if np.any(mask_b):
            cols_b = cols[mask_b]
            q_b = q_arr[mask_b]
            other = ll_other_arr[q_b]
            delta = ll_this_arr[q_b] - other

            llA[cols_b] += other
            llT[cols_b] += other
            llC[cols_b] += other
            llG[cols_b] += other

            idx = base_to_idx[base_arr[mask_b]]
            if np.any(idx == 0):
                llA[cols_b[idx == 0]] += delta[idx == 0]
            if np.any(idx == 1):
                llT[cols_b[idx == 1]] += delta[idx == 1]
            if np.any(idx == 2):
                llC[cols_b[idx == 2]] += delta[idx == 2]
            if np.any(idx == 3):
                llG[cols_b[idx == 3]] += delta[idx == 3]

    ll_sums = np.stack([llA, llT, llC, llG], axis=0)

    full_ll = logsumexp(ll_sums, axis=0)

    prob_max = np.asarray(np.exp(np.amax(ll_sums, axis=0) - full_ll)).ravel()
    nuc_max = np.asarray(np.argmax(ll_sums, axis=0)).ravel()

    nuc_arr = np.array(nucleotides, dtype='<U1')
    seq_arr = np.where(prob_max > 0.3, nuc_arr[nuc_max], 'N')
    master_read['seq'] = ''.join(seq_arr.tolist())
    master_read['phred'] = np.nan_to_num(np.rint(-10*np.log10(1-prob_max+1e-13)))

    if len(reverse_read1) == 0:
        return (False, ':'.join([gene,cell,umi]))
    v, c = np.unique(reverse_read1, return_counts=True)
    m = c.argmax()

    master_read['SN'] = read.reference_name
    master_read['is_reverse'] = v[m]
    master_read['ref_intervals'] = interval(intervals_extract(np.sort(ref_pos_set_array)))
    master_read['skipped_intervals'] = interval(list(set(master_read['skipped_intervals'])))
    master_read['del_intervals'] =  ~(master_read['ref_intervals'] | master_read['skipped_intervals'])
    master_read['NR'] = nreads
    master_read['IR'] = np.sum(intronic_list)
    master_read['ER'] = np.sum(exonic_list)
    master_read['cell'] = cell
    master_read['gene'] = gene
    master_read['umi'] = umi
    return (True, convert_to_sam(master_read, UMI_tag))

def get_compatible_isoforms_stitcher(mol_list, isoform_dict_json,refskip_dict_json, h):
    def interval_total_len(x):
        return sum((t[2] - t[1] + 1) for t in P.to_data(x))

    isoform_dict = P.IntervalDict()
    for i,s in isoform_dict_json.items():
        isoform_dict[P.from_string(i, conv=int)] = set(s.split(','))
    refskip_dict = P.IntervalDict()
    for i,s in refskip_dict_json.items():
        refskip_dict[P.from_string(i, conv=int)] = set(s.split(','))
    
    compatible_isoforms_trie = dict()
    new_mol_list = []
    for success, m in mol_list:
        if not success:
            if type(m) is str:
                new_mol_list.append((success,m))
            else:
                new_mol_list.append((success,m.to_string()))
            continue
        mol = pysam.AlignedRead.fromstring(m,h)
        i = interval(intervals_extract(mol.get_reference_positions()))
        refskip_cigar = [t[0] for t in mol.cigartuples if t[1] > 0 and t[0] in [2,3]]
        blocks = mol.get_blocks()
        j = []
        for n in range(len(blocks)-1):
            if refskip_cigar[n] == 3:
                j.append((blocks[n][1],blocks[n+1][0]))
        j = interval(j)
        set_list = [s for k,s in isoform_dict.get(i, default={'intronic'}).items() if interval_total_len(k) > 4]
        set_refskip_list = [s for k,s in refskip_dict.get(j, default={'intronic'}).items() if interval_total_len(k) > 4]
        if {'intronic'} in set_list:
            if len(set_list) > 1:
                del set_list[set_list.index({'intronic'})]
        if {'intronic'} in set_refskip_list:
            if len(set_refskip_list) > 1:
                del set_refskip_list[set_refskip_list.index({'intronic'})]
        try:
            if len(set_refskip_list) > 0:
                mol.set_tag('CT',','.join(list(set.intersection(*set_list).intersection(*set_refskip_list))))
            else:
                mol.set_tag('CT',','.join(list(set.intersection(*set_list))))
            new_mol_list.append((success,mol.to_string()))
        except:
            continue
    return new_mol_list

def assemble_reads(bamfile,gene_to_stitch, cell_set, isoform_dict_json,refskip_dict_json,single_end, cell_tag, UMI_tag, q):
    read_dict = {}
    bam = pysam.AlignmentFile(bamfile, 'rb')
    gene_of_interest = gene_to_stitch['gene_id']
    for read in bam.fetch(gene_to_stitch['seqid'], gene_to_stitch['start'], gene_to_stitch['end']):
        if read.has_tag('CB'):
            cell = read.get_tag('CB')
        else:
            continue
        if cell_set is not None:
            if cell not in cell_set:
                continue
        umi = read.get_tag(UMI_tag)
        if umi == '':
            continue
        else:
            if read.has_tag('GE'):
                gene_exon = read.get_tag('GE')
            else:
                gene_exon = 'Unassigned'
            if read.has_tag('GI'):
                gene_intron = read.get_tag('GI')
            else:
                gene_intron = 'Unassigned'
            # if it maps to the intron or exon of a gene
            if gene_intron != 'Unassigned' or gene_exon != 'Unassigned':
                # if it is a junction read
                if gene_intron == gene_exon:
                    gene = gene_intron
                    # if it's an only intronic read
                elif gene_intron != 'Unassigned' and gene_exon == 'Unassigned':
                    gene = gene_intron
                    # if it's an only exonic read
                elif gene_exon != 'Unassigned' and gene_intron == 'Unassigned':
                    gene = gene_exon
                    # if the exon and intron gene tag contradict each other
                else:
                    continue
            else:
                continue
        if single_end:
            if gene == gene_of_interest and not read.is_unmapped:
                node = (cell, gene, umi)
                if node in read_dict:
                    read_dict[node].append(read)
                else:
                    read_dict[node] = [read]
        else:
            if read.is_paired and not read.is_unmapped and not read.mate_is_unmapped and gene == gene_of_interest and read.is_proper_pair:
                node = (cell, gene, umi)
                if node in read_dict:
                    read_dict[node].append(read)
                else:
                    read_dict[node] = [read]
    mol_list = []
    mol_append = mol_list.append
    for node, mol in read_dict.items():
        n_read1 = np.sum([(r.is_read1)&(r.get_tag(UMI_tag) != '') for r in mol])
        if n_read1 > 0:
            mol_append(stitch_reads(mol, single_end, node[0], node[1], node[2], UMI_tag))
    del read_dict
    if isoform_dict_json is not None:
        mol_list = get_compatible_isoforms_stitcher(mol_list, isoform_dict_json,refskip_dict_json, bam.header)
    if len(mol_list) == 0:
        return gene_of_interest
    if len(mol_list) > 50000:
        for m_list in chunks(mol_list, 50000):
            q.put((True, m_list))
    else:
        q.put((True, mol_list))
    return gene_of_interest


def make_POS_and_CIGAR(stitched_m):
    CIGAR = ''
    conflict = False
    interval_list = []
    ref_and_skip_intersect = stitched_m['ref_intervals'] & stitched_m['skipped_intervals']
    nreads_conflict = 0
    if not ref_and_skip_intersect.empty:
        conflict = True
        nreads_conflict = len(list(P.iterate(ref_and_skip_intersect, step=1))) 
        stitched_m['skipped_intervals'] = stitched_m['skipped_intervals'] - ref_and_skip_intersect
        interval_list = [i for t in P.to_data(ref_and_skip_intersect) for i in t[1:-1]]
    ref_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['ref_intervals'])]
    if stitched_m['skipped_intervals'].empty:
        skipped_tuples = []
    else:
        skipped_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['skipped_intervals'])]
    if stitched_m['del_intervals'].empty:
        del_tuples = []
    else:
        del_tuples = [(i[1] if i[0] else i[1]+1, i[2] if i[3] else i[2]-1) for i in P.to_data(stitched_m['del_intervals'])[1:-1]]
    POS = ref_tuples[0][0] + 1
    tuple_dict = {'M': ref_tuples, 'N': skipped_tuples, 'D': del_tuples}
    while sum(len(t) for t in tuple_dict.values()) > 0:
        pos_dict = {k:v[0][0] for k,v in tuple_dict.items() if len(v) > 0}
        c = min(pos_dict, key=pos_dict.get)
        n_bases = np.int_(tuple_dict[c[0]][0][1]-tuple_dict[c[0]][0][0])+1
        if n_bases == 0:
            del tuple_dict[c[0]][0]
            continue
        CIGAR += '{}{}'.format(n_bases,c[0])
        del tuple_dict[c[0]][0]
    return POS, CIGAR, conflict, nreads_conflict, interval_list

def convert_to_sam(stitched_m, UMI_tag):
    sam_dict = {}
    POS, CIGAR, conflict, nreads_conflict, interval_list = make_POS_and_CIGAR(stitched_m)
    sam_dict['QNAME'] = '{}:{}:{}'.format(stitched_m['cell'],stitched_m['gene'],stitched_m['umi'])
    sam_dict['FLAG'] = str(16*stitched_m['is_reverse'])
    sam_dict['RNAME'] = stitched_m['SN']
    sam_dict['POS'] = str(POS)
    sam_dict['MAPQ'] = str(255)
    sam_dict['CIGAR'] = CIGAR
    sam_dict['RNEXT'] = '*'
    sam_dict['PNEXT'] = str(0)
    sam_dict['TLEN'] = str(0)
    sam_dict['SEQ'] = stitched_m['seq']
    sam_dict['QUAL'] = "".join([chr(int(p)) for p in np.clip(stitched_m['phred'],0,126-33)+33])
    sam_dict['NR'] = 'NR:i:{}'.format(stitched_m['NR'])
    sam_dict['ER'] = 'ER:i:{}'.format(stitched_m['ER'])
    sam_dict['IR'] = 'IR:i:{}'.format(stitched_m['IR'])
    sam_dict['BC'] = 'BC:Z:{}'.format(stitched_m['cell'])
    sam_dict['XT'] = 'XT:Z:{}'.format(stitched_m['gene'])
    sam_dict[UMI_tag] = '{}:Z:{}'.format(UMI_tag, stitched_m['umi'])
    #sam_dict['EL'] = 'EL:B:I,{}'.format(','.join([str(e) for e in stitched_m['ends']]))
    if conflict:
        sam_dict['NC'] = 'NC:i:{}'.format(nreads_conflict)
        sam_dict['IL'] = 'IL:B:I,{}'.format(','.join([str(e) for e in interval_list]))
    return '\t'.join(list(sam_dict.values()))

def yield_reads(read_dict):
    for cell in read_dict:
        for gene in read_dict[cell]:
            #print('\t', gene)
            for umi in read_dict[cell][gene]:
                #print('\t\t', umi)
                yield read_dict[cell][gene][umi], None, cell, gene, umi


def create_write_function(filename, bamfile, version):
    bam = pysam.AlignmentFile(bamfile, 'rb')
    header = bam.header
    
    def write_sam_file(q):
        error_file = open('{}_error.log'.format(os.path.splitext(filename)[0]), 'w')
        stitcher_bam = pysam.AlignmentFile(filename,'wb',header={'HD':header['HD'], 'SQ':header['SQ'], 'PG': [{'ID': 'stitcher.py','VN': '{}'.format(version)}]})
        while True:
            good, mol_list = q.get()
            if good is None: break
            if good:
                g = ''
                for success, mol in mol_list:
                    if success:
                        read = pysam.AlignedRead.fromstring(mol,header)
                        if g == '':
                            g = read.get_tag('XT')
                        stitcher_bam.write(read)
                    else:
                        error_file.write(mol+'\n')
                if g != '':
                    error_file.write('Gene:{}\n'.format(g))
            q.task_done()
        q.task_done()
        error_file.close()
        stitcher_bam.close()
        return None
    return write_sam_file

def extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)
    
def construct_stitched_molecules(infile, outfile,gtffile,isoformfile, junctionfile, cells, gene_file, contig, threads, single_end, cell_tag, UMI_tag, gene_identifier, skip_iso, q, version):
    if cells is not None:
        cell_set = set([line.rstrip() for line in open(cells)])
    else:
        cell_set = None
    print('Reading gene info from {}'.format(gtffile))
    gene_list = []
    if gene_identifier == 'gene_id':
        n = 1
    elif gene_identifier == 'gene_name':
        n = 5
    else:
        n = 1
    with open(gtffile, 'r') as f:
        for line in f:
            l = line.split('\t')
            if len(l) < 8:
                continue
            if l[2] == 'gene':
                if contig is not None:
                    if l[0] == contig:
                        try:
                            gene_list.append({'gene_id': l[8].split(' ')[n].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
                        except:
                            gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
                    else:
                        continue
                else:
                    try:
                        gene_list.append({'gene_id': l[8].split(' ')[n].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
                    except:
                        gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';\n'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4])})
    gene_dict = {g['gene_id']: g for g in gene_list}
    
    if gene_file is not None and gene_file != 'None':
        gene_set = set([line.rstrip() for line in open(gene_file)])
        gene_dict = {k:v for k,v in gene_dict.items() if k in gene_set}
    bam = pysam.AlignmentFile(infile, 'rb')
    contig_set = set([d['SN'] for d in bam.header['SQ']])
    prev_l = len(gene_dict)
    gene_dict = {k:v for k,v in gene_dict.items() if v['seqid'] in contig_set}
    new_l = len(gene_dict)
    diff_l = prev_l - new_l
    print(prev_l, new_l)
    if diff_l > 0:
        warnings.warn('Warning: removed {diff_l} genes with contig not present in bam file'.format(diff_l=diff_l))
    bam.close()
    if skip_iso:
        print('Skipping isoform info')
        params = Parallel(n_jobs=threads, verbose = 3, backend='loky')(delayed(assemble_reads)(infile, gene, cell_set,None, None, single_end, cell_tag, UMI_tag, q) for g,gene in gene_dict.items())
    else:    
        print('Reading isoform info from {}'.format(isoformfile))
        with open(isoformfile) as json_file:
            isoform_unique_intervals = json.load(json_file)
        with open(junctionfile) as json_file:
            refskip_unique_intervals = json.load(json_file)
        params = Parallel(
            n_jobs=threads, verbose = 3, backend='loky'
            )(
                delayed(assemble_reads)(
                    infile, gene, cell_set,isoform_unique_intervals[g],refskip_unique_intervals[g],single_end, cell_tag, UMI_tag, q
                    ) for g,gene in gene_dict.items() if g in isoform_unique_intervals and g in refskip_unique_intervals
                    )


    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stitch together molecules from reads sharing the same UMI')
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', metavar='output', type=str, help='Output .bam file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('-iso','--isoform',metavar='iso', type=str, help='json file with isoform information')
    parser.add_argument('-jun','--junction', metavar='jun', type=str, help='json file with exon-exon structure')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--single-end', action='store_true', help='Activate flag if data is single-end')
    parser.add_argument('--skip-iso', action='store_true', help='Skip isoform calling')
    parser.add_argument('--UMI-tag', type=str, default='UB', help='UMI tag to group reads')
    parser.add_argument('--cell-tag', type=str, default='BC', help='cell baroced tag to group reads')
    parser.add_argument('--cells', default=None, metavar='cells', type=str, help='List of cell barcodes to stitch molecules')
    parser.add_argument('--genes', default=None, metavar='genes', type=str, help='List of gene,  one per line.')
    parser.add_argument('--contig', default=None, metavar='contig', type=str, help='Restrict stitching to contig')
    parser.add_argument('--gene-identifier', default='gene_id', metavar='gene_identifier', type=str, help='Gene identifier (gene_id or gene_name)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()
    infile = args.input
    if infile is None:
        raise Exception('No input file provided.')
    outfile = args.output
    if outfile is None:
        raise Exception('No output file provided.')
    gtffile = args.gtf  
    if gtffile is None:
        raise Exception('No gtf file provided.')
    skip_iso = args.skip_iso
    if not skip_iso:
        isoformfile = args.isoform
        junctionfile = args.junction
    else:
        isoformfile = ''
        junctionfile = ''
    threads = int(args.threads)
    cells = args.cells
    gene_file = args.genes
    contig = args.contig
    single_end = args.single_end
    UMI_tag = args.UMI_tag
    cell_tag = args.cell_tag
    gene_identifier = args.gene_identifier
    m = Manager()
    q = m.JoinableQueue()
    p = Process(target=create_write_function(filename=outfile, bamfile=infile, version=__version__), args=(q,))
    p.start()
    
    print('Stitching reads for {}'.format(infile))
    
    start = time.time()
    construct_stitched_molecules(infile, outfile, gtffile, isoformfile,junctionfile, cells, gene_file, contig, threads,single_end,cell_tag, UMI_tag,gene_identifier, skip_iso, q, __version__)
    q.put((None,None))
    p.join()
    end = time.time()
    
    print('Finished writing stitched molecules from {} to {}, took {}'.format(infile, outfile, get_time_formatted(end-start)))
