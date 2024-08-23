options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("./src/collate_fastq_reads.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("-c", "--uncollapsed_raw_counts", default="raw_counts_uncollapsed.csv",
                    help="path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("--sequencing_index_cols", default= "index_1,index_2", 
                    help = "Sequencing columns in the sample meta")
parser$add_argument("--id_cols", default= "pcr_plate,pcr_well", 
                    help = "Columns that identify a unique PCR well")
parser$add_argument("--reverse_index2", action="store_true", default=FALSE, 
                    help= "Reverse complement of index 2 for NovaSeq and NextSeq")
parser$add_argument("--barcode_col", default= "forward_read_cl_barcode", 
                    help= "Name of the column in uncollapsed_raw_counts that contains the barcode sequences.")
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
  # Read in files and parse vector arguments
  uncollapsed_raw_counts= data.table::fread(expected_file_path, header= T, sep= ',', data.table= F)
  sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',', data.table= F)
  
  # Parse vector inputs
  sequencing_index_cols= unlist(strsplit(args$sequencing_index_cols, ","))
  id_cols= unlist(strsplit(args$id_cols, ","))
  
  # Validation: Check that sequencing_index_cols are from sample meta column names
  if(!all(sequencing_index_cols %in% colnames(sample_meta))) {
    stop(paste('The following sequencing_index_cols were not found in the sample meta: ',
               sequencing_index_cols[!sequencing_index_cols %in% colnames(sample_meta)]))
  }
  
  # Validation: Check that id_cols are from sample meta column names
  if(!all(id_cols %in% colnames(sample_meta))) {
    stop(paste('The following id_cols were not found in the sample meta: ',
               id_cols[!id_cols %in% colnames(sample_meta)]))
  }
  
  print("Collating fastq reads ...")
  raw_counts= collate_fastq_reads(uncollapsed_raw_counts, sample_meta, 
                                  sequencing_index_cols= sequencing_index_cols,
                                  id_cols= id_cols,
                                  reverse_index2= args$reverse_index2,
                                  barcode_col= args$barcode_col)
  
  # Validation: Basic file size check
  if(nrow(raw_counts) == 0) {
    stop('ERROR: Empty file generated. No rows in raw_counts output.')
  } 
  
  rc_out_file= paste(args$out, 'raw_counts.csv', sep='/')
  print(paste("Writing raw_counts.csv to ", rc_out_file))
  write.csv(raw_counts, rc_out_file, row.names= FALSE, quote= FALSE)
} else {
  print("Uncollapsed raw counts file not detected. Proceeding with generating filtered counts file.")
}
