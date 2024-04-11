library(argparse)
library(magrittr)
library(tidyverse)
source("./src/collate_fastq_reads.R")

parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--uncollapsed_raw_counts", default="raw_counts_uncollapsed.csv",
                    help="path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

#Save arguments
if (args$out == ""){
  args$out = args$wkdir
}

if (file.exists(args$uncollapsed_raw_counts)) {
  sample_meta = read.csv(args$sample_meta)
  uncollapsed_raw_counts = read.csv(args$uncollapsed_raw_counts)
  print("Collating fastq reads")
  raw_counts <- collate_fastq_reads(sample_meta,uncollapsed_raw_counts) 
  
  rc_out_file = paste(
    args$out,
    'raw_counts.csv',
    sep='/'
  )
  print(paste("Writing to file: ", rc_out_file))
  write.csv(raw_counts, rc_out_file, row.names=F, quote=F)
} else {
  print("Uncollapsed raw counts file not detected. Proceeding with generating filtered counts file.")
}
