library(argparse)
library(magrittr)
library(tidyverse)
source("./src/collate_fastq_reads.R")

# Parser ----
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("-c", "--uncollapsed_raw_counts", default="raw_counts_uncollapsed.csv",
                    help="path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("--seq_cols", default= "IndexBarcode1,IndexBarcode2", help = "Sequencing columns in the sample meta")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == "") {
  args$out = args$wkdir
}

# Run collate_fastq_reads if uncollapsed file exists ----
expected_file_path <- paste(args$out, "raw_counts_uncollapsed.csv", sep='/')

if(file.exists(expected_file_path)) {
  sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',', data.table= F)
  uncollapsed_raw_counts= data.table::fread(expected_file_path, header= T, sep= ',', data.table= F)
  seq_cols= unlist(strsplit(args$seq_cols, ","))
  
  # QC: Check if seq_cols is composed of sample meta column names
  if (!all(seq_cols %in% colnames(sample_meta))) {
    stop(paste("Colnames not found in sample_meta, check metadata or --seq_cols argument:", args$seq_cols))
  }
  
  print("Collating fastq reads")
  raw_counts= collate_fastq_reads(uncollapsed_raw_counts, sample_meta, seq_cols)
  
  # QC: Basic file size check
  if(nrow(raw_counts) == 0) {
    stop('ERROR: Empty file generated. No rows in raw_counts output.')
  } 
  
  rc_out_file = paste(args$out, 'raw_counts.csv', sep='/')
  print(paste("Writing to file: ", rc_out_file))
  write.csv(raw_counts, rc_out_file, row.names=F, quote=F)
} else {
  print("Uncollapsed raw counts file not detected. Proceeding with generating filtered counts file.")
}
