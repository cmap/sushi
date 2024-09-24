options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
source("./src/collate_fastq_reads.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument('--raw_counts_uncollapsed', default="raw_counts_uncollapsed.csv",
                    help="path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("--sequencing_index_cols", default= "index_1,index_2", 
                    help = "Sequencing columns in the sample meta")
parser$add_argument("--id_cols", default= "pcr_plate,pcr_well", 
                    help = "Columns that identify a unique PCR well")
parser$add_argument("--reverse_index2", type="logical", default=FALSE,
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

# Read in sample meta and parse argument strings ----
# Read in files and parse vector arguments
sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',')

# Parse vector inputs
sequencing_index_cols= unlist(strsplit(args$sequencing_index_cols, ","))
id_cols= unlist(strsplit(args$id_cols, ","))

# Validation: Check that sequencing_index_cols present in the sample meta ----
if(!validate_columns_exist(sequencing_index_cols, sample_meta)) {
  print('The following sequencing_index_cols are not present in the sample meta.')
  stop('One or more sequencing_index_cols is NOT present in the sample meta.')
}

# Validation: Check that id_cols are present in the sample meta ----
if(!validate_columns_exist(id_cols, sample_meta)) {
  stop('One or more id_cols is NOT present in the sample meta.')
}

# Run collate_fastq_reads on chunks of raw_counts_uncollapsed.csv ----
summed_reads= process_in_chunks(large_file_path= args$raw_counts_uncollapsed, 
                                chunk_size= 10^6, 
                                action= collate_fastq_reads,
                                sample_meta= sample_meta, 
                                sequencing_index_cols= sequencing_index_cols,
                                id_cols= id_cols,
                                reverse_index2= args$reverse_index2,
                                barcode_col= args$barcode_col)

# Sum up the read across the chunks afterwards!
summed_reads= summed_reads[, .(n= sum(n)), by= c(id_cols, args$barcode_col)]

# Split reads by either known or unknown ----
# Reads are separated by whether or not the barcode exists in the PRISM library
# Read in metadata to get list of all known barcodes
cell_line_meta= data.table::fread(args$cell_line_meta, header= TRUE, sep= ',')
CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')

# Call function to separate barcodes
split_reads= extract_known_barcodes(summed_reads, unique(c(cell_line_meta$Sequence, CB_meta$Sequence)),
                                    barcode_col= args$barcode)

# Validation: Basic file size check ----
if(nrow(split_reads$mapped_reads) == 0) {
  stop('ERROR: Empty file generated. No rows in raw_counts output.')
} 

# Write out files ----
out_file= paste(args$out, 'unknown_reads.csv', sep='/')
print(paste("Writing unknown_reads.csv to ", out_file))
write.csv(split_reads$unknown_reads, out_file, row.names= FALSE, quote= FALSE)

out_file= paste(args$out, 'known_reads.csv', sep='/')
print(paste("Writing known_reads.csv to ", out_file))
write.csv(split_reads$known_reads, out_file, row.names= FALSE, quote= FALSE)
