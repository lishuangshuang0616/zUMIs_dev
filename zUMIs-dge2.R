#!/usr/bin/env Rscript
library(methods)
library(data.table)
library(yaml)
library(ggplot2)
suppressPackageStartupMessages(library(Rsamtools))

##########################
myYaml <- commandArgs(trailingOnly = T)

opt   <-read_yaml(myYaml)
setwd(opt$out_dir)
#try(unixtools::set.tempdir(opt$out_dir))
source(paste0(opt$zUMIs_directory,"/runfeatureCountFUN.R"))
source(paste0(opt$zUMIs_directory,"/misc/featureCounts.R"))
source(paste0(opt$zUMIs_directory,"/UMIstuffFUN.R"))
source(paste0(opt$zUMIs_directory,"/barcodeIDFUN.R"))
options(datatable.fread.input.cmd.message=FALSE)
print(Sys.time())

samtoolsexc <- opt$samtools_exec
data.table::setDTthreads(threads=1)

#Check the version of Rsubread
#checkRsubreadVersion()
fcounts_clib <- paste0(opt$zUMIs_directory,"/misc/fcountsLib2")

opt <- fixMissingOptions(opt)
#######################################################################
########################## double check for non-UMI method
UMIcheck <- check_nonUMIcollapse(opt$sequence_files)
if(UMIcheck == "nonUMI"){
  opt$counting_opts$Ham_Dist <- 0
}
#is the data Smart-seq3?
smart3_flag <- ifelse(any(grepl(pattern = "ATTGCGCAATG",x = unlist(opt$sequence_files))), TRUE, FALSE)

#######################################################################
##### Barcode handling & chunking

#read file with barcodecounts

#check if binning of adjacent barcodes should be run
if(opt$barcodes$BarcodeBinning > 0){
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes_binned.txt"))
}else{
  bccount <- fread(paste0(opt$out_dir,"/zUMIs_output/",opt$project,"kept_barcodes.txt"))
}
bccount<-splitRG(bccount=bccount, mem= opt$mem_limit, hamdist = opt$counting_opts$Ham_Dist)

##############################################################
##### featureCounts

abamfile<-paste0(opt$out_dir,"/",opt$project,".filtered.tagged.Aligned.out.bam")
outbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.bam")

## gene annotation
saf<-.makeSAF(gtf = paste0(opt$out_dir,"/",opt$project,".final_annot.gtf"),
              extension_var = opt$reference$exon_extension,
              exon_extension = opt$reference$extension_length,
              buffer_length = (opt$reference$extension_length / 2),
              scaff_length = opt$reference$scaffold_length_min,
              multi_overlap_var = opt$counting_opts$multi_overlap,
              samtoolsexc = samtoolsexc)
try(gene_name_mapping <- .get_gene_names(gtf = paste0(opt$out_dir,"/",opt$project,".final_annot.gtf")), silent = TRUE)
try(data.table::fwrite(gene_name_mapping, file = paste0(opt$out_dir,"/zUMIs_output/expression/",opt$project,".gene_names.txt"), sep ="\t", quote = FALSE), silent = TRUE)
##

if(opt$counting_opts$strand == 1){
  #split bam in UMU ends and internal
  print("Preparing Smart-seq3 data for stranded gene assignment...")
  print(Sys.time())
  tmp_bams <- split_bam(bam = abamfile, cpu = opt$num_threads, samtoolsexc=samtoolsexc)

  #assign features with appropriate strand
  fnex_int<-.runFeatureCount(tmp_bams[1], saf=saf$exons, strand=0, type="ex", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap, fraction_overlap = opt$counting_opts$fraction_overlap)
  fnex_umi<-.runFeatureCount(tmp_bams[2], saf=saf$exons, strand=1, type="ex", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap, fraction_overlap = opt$counting_opts$fraction_overlap)
  ffiles_int <- paste0(fnex_int,".tmp")
  ffiles_umi <- paste0(fnex_umi,".tmp")

  if(opt$counting_opts$introns){
    fnin_int<-.runFeatureCount(ffiles_int, saf=saf$introns, strand=0, type="in", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap, fraction_overlap = opt$counting_opts$fraction_overlap)
    fnin_umi<-.runFeatureCount(ffiles_umi, saf=saf$introns, strand=1, type="in", primaryOnly = opt$counting_opts$primaryHit, cpu = opt$num_threads, mem = opt$mem_limit, fcounts_clib = fcounts_clib, multi_overlap_var = opt$counting_opts$multi_overlap, fraction_overlap = opt$counting_opts$fraction_overlap)
    ffiles_int <- paste0(fnin_int,".tmp")
    ffiles_umi <- paste0(fnin_umi,".tmp")
  }
  join_bam_cmd <- paste(samtoolsexc, "cat -o", outbamfile, ffiles_int, ffiles_umi)
  system(join_bam_cmd)
  system(paste0("rm ",tmp_bams[1],"* ",tmp_bams[2],"*"))
}else{
  fnex<-.runFeatureCount(abamfile,
                         saf=saf$exons,
                         strand=opt$counting_opts$strand,
                         type="ex",
                         primaryOnly = opt$counting_opts$primaryHit,
                         cpu = opt$num_threads,
                         mem = opt$mem_limit,
                         fcounts_clib = fcounts_clib,
                         multi_overlap_var = opt$counting_opts$multi_overlap,
                         fraction_overlap = opt$counting_opts$fraction_overlap)
  ffiles<-paste0(fnex,".tmp")

  if(opt$counting_opts$introns){
    fnin  <-.runFeatureCount(ffiles,
                             saf=saf$introns,
                             strand=opt$counting_opts$strand,
                             type="in",
                             primaryOnly = opt$counting_opts$primaryHit,
                             cpu = opt$num_threads,
                             mem = opt$mem_limit,
                             fcounts_clib = fcounts_clib,
                             multi_overlap_var = opt$counting_opts$multi_overlap,
                             fraction_overlap = opt$counting_opts$fraction_overlap)
    system(paste0("rm ",fnex,".tmp"))
    ffiles<-paste0(fnin,".tmp")
  }

  system(paste("mv",ffiles,outbamfile))
}

if(is.null(opt$mem_limit)){
  mempercpu <- max(round(100/opt$num_threads,0),1)
}else{
  mempercpu <- max(round(opt$mem_limit/opt$num_threads,0),1)
}


# === PYTHON HANDOVER ===
# The following R logic for sorting, UMI collapsing, and counting is disabled.
# It is now handled by dge_analysis.py for better performance.
print("R script finished FeatureCounts. Handing over to Python for UMI processing.")
print(Sys.time())
q()

# OLD CODE COMMENTED OUT FOR REFERENCE
# if(opt$counting_opts$Ham_Dist == 0){
#   sortbamfile <-paste0(opt$out_dir,"/",opt$project,".filtered.Aligned.GeneTagged.sorted.bam")
# ... (rest of the file)

