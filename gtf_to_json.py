import portion as P
import itertools
from joblib import delayed, Parallel
import gffutils
import json
import argparse
import os
import sys

def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def interval(t):
    return P.from_data([(True,i[0],i[1], True) for i in t])

def create_interval_dict_linear_time(gene, isoform_interval_dict):
    d = P.IntervalDict()
    events = []
    for transcript, inter in isoform_interval_dict.items():
        for lower_closed, lower, upper, upper_closed in P.to_data(inter):
            start = int(lower) + (0 if lower_closed else 1)
            end = int(upper) - (0 if upper_closed else 1)
            if end < start:
                continue
            events.append((start, transcript, 1))
            events.append((end + 1, transcript, -1))

    if not events:
        return gene, d

    events.sort(key=lambda x: x[0])
    active = set()
    last_pos = events[0][0]
    segs_by_set = {}

    i = 0
    n_events = len(events)
    while i < n_events:
        pos = events[i][0]
        if pos > last_pos and active:
            key = tuple(sorted(active))
            if key not in segs_by_set:
                segs_by_set[key] = []
            segs_by_set[key].append((last_pos, pos - 1))
        
        while i < n_events and events[i][0] == pos:
            _, tr, delta = events[i]
            if delta > 0:
                active.add(tr)
            else:
                active.discard(tr)
            i += 1
        last_pos = pos

    for tr_set, segs in segs_by_set.items():
        # Merging segments is only needed if they are fragmented, 
        # but our sweep line produces contiguous chunks for identical active sets unless interrupted.
        # However, to be safe and match original logic:
        segs.sort()
        merged = []
        for s, e in segs:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        d[interval(merged)] = set(tr_set)
    return gene, d

def process_gene_data(gene_id, transcripts_data):
    """
    Process a single gene's transcript data to calculate unique intervals.
    transcripts_data: dict {transcript_id: [(start, end), ... sorted exons]}
    """
    isoform_interval_dict = {}
    isoform_refskip_dict = {}

    for t_id, exons in transcripts_data.items():
        # Build Exon Intervals
        # Using P.closed(start, end) for each exon
        # Performance: P.empty() | P.closed(...)
        curr_intervals = P.empty()
        for start, end in exons:
            curr_intervals |= P.closed(start, end)
        isoform_interval_dict[t_id] = curr_intervals

        # Build Refskip Intervals (Introns)
        # Genomic gap between exons
        curr_refskip = P.empty()
        # Exons are sorted by start (genomic coordinates)
        # Regardless of strand, the intron is the gap between exon[i].end and exon[i+1].start
        for i in range(len(exons) - 1):
            # Using P.closed to include boundaries as per original script logic (implied)
            # or perhaps P.open/closed?
            # Original code: P.closed(exons[i].end, exons[i+1].start)
            # Assuming boundaries are inclusive 1-based gffutils coords.
            # Intron starts after exon end and ends before next exon start.
            # If we strictly follow original code's positive strand logic:
            curr_refskip |= P.closed(exons[i][1], exons[i+1][0])
        
        isoform_refskip_dict[t_id] = curr_refskip

    # Calculate Unique Intervals
    _, unique_intervals_dict = create_interval_dict_linear_time(gene_id, isoform_interval_dict)
    _, unique_refskips_dict = create_interval_dict_linear_time(gene_id, isoform_refskip_dict)

    # Serialize for JSON
    # Format: { interval_string : comma_joined_transcripts }
    ui_dump = {P.to_string(k): ','.join(sorted(v)) for k, v in unique_intervals_dict.items()}
    ur_dump = {P.to_string(k): ','.join(sorted(v)) for k, v in unique_refskips_dict.items()}

    return gene_id, ui_dump, ur_dump

def fetch_gene_data(db):
    """Generator that yields gene data for processing."""
    for gene in db.features_of_type('gene'):
        g_id = gene['gene_id'][0]
        transcripts_data = {}
        
        # db.children allows efficient retrieval
        # Retrieve all exons for the gene grouped by transcript
        # To avoid multiple queries, we can query all transcripts first
        
        for transcript in db.children(gene, featuretype='transcript'):
            t_id = transcript['transcript_id'][0]
            # Get exons for this transcript
            exons = []
            for exon in db.children(transcript, featuretype='exon'):
                exons.append((exon.start, exon.end))
            
            if exons:
                exons.sort() # Ensure sorted by genomic start
                transcripts_data[t_id] = exons
        
        if transcripts_data:
            yield g_id, transcripts_data

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Write json file used for stitcher.py from a gtf file')
    parser.add_argument('-g','--gtf',metavar='gtf', type=str, help='Input gtf file')
    parser.add_argument('-d','--db', metavar='db', type=str, help='Intermediary database (db) file')
    parser.add_argument('-ji','--json_intervals', metavar='json', type=str, help='Output json file for coverage')
    parser.add_argument('-jr','--json_refskip', metavar='json', type=str, help='Output json file for refskip')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--force-db', action='store_true', help='Recreate database even if it exists')
    args = parser.parse_args()
    
    gtffile = args.gtf
    dbfile = args.db
    jsonfile_1 = args.json_intervals
    jsonfile_2 = args.json_refskip
    threads = int(args.threads)

    print('Creating/Loading gtf database...')
    # Use gffutils to create or load DB
    # Disable check_sq for speed if not strictly needed
    if dbfile and os.path.exists(dbfile) and not args.force_db:
        db = gffutils.FeatureDB(dbfile, keep_order=True)
    else:
        db = gffutils.create_db(
            gtffile,
            dbfile,
            force=True,
            keep_order=True,
            merge_strategy='merge',
            sort_attribute_values=True,
            disable_infer_genes=True, # Speed up if genes are explicit
            disable_infer_transcripts=True # Speed up
        )

    print('Loading gene data into memory to avoid SQLite threading issues...')
    # Materialize the generator to a list in the main thread
    # This ensures all DB access happens here, before joblib threads start
    gene_data_list = list(fetch_gene_data(db))
    
    print(f'Loaded {len(gene_data_list)} genes. Starting parallel processing...')
    
    # Use joblib to process genes in parallel
    results = Parallel(n_jobs=threads, verbose=3, backend='loky')(
        delayed(process_gene_data)(g_id, t_data) 
        for g_id, t_data in gene_data_list
    )

    # Aggregate results
    isoform_unique_intervals_dump = {}
    isoform_unique_refskip_dump = {}

    for gid, ui, ur in results:
        if ui:
            isoform_unique_intervals_dump[gid] = ui
        if ur:
            isoform_unique_refskip_dump[gid] = ur

    print('Writing unique isoform intervals to json file {}'.format(jsonfile_1))
    with open(jsonfile_1, 'w') as fp:
        json.dump(isoform_unique_intervals_dump, fp)

    print('Writing unique isoform refskip to json file {}'.format(jsonfile_2))
    with open(jsonfile_2, 'w') as fp:
        json.dump(isoform_unique_refskip_dump, fp)

