options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
source("./src/collate_fastq_reads.R")
source("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument('--raw_counts_uncollapsed', default= "raw_counts_uncollapsed.csv",
                    help= "path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("--cell_line_meta", default="cell_line_meta.csv", help= "Cell line metadata")
parser$add_argument("--CB_meta", default= "CB_meta.csv", help= "Control Barcode metadata")
parser$add_argument('--sequencing_index_cols', default= 'index_1,index_2', 
                    help= 'List of sequencing columns in the sample meta.')
parser$add_argument("--id_cols", default= "pcr_plate,pcr_well", 
                    help = "Columns that identify a unique PCR well")
parser$add_argument("--reverse_index2", type="logical", default=FALSE,
                    help= "Reverse complement of index 2 for NovaSeq and NextSeq")
parser$add_argument("--barcode_col", default= "forward_read_cl_barcode", 
                    help= "Name of the column in uncollapsed_raw_counts that contains the barcode sequences.")
parser$add_argument('--low_abundance_threshold', default= 20, 
                    help= 'For unknown barcodes, counts below this threshold will be marked as an unknown barcode.')
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == "") {
  args$out = args$wkdir
}

# Read in metadata files as data.table objects ----
sample_meta= data.table::fread(args$sample_meta, header= TRUE, sep= ',')
cell_line_meta= data.table::fread(args$cell_line_meta, header= TRUE, sep= ',')
CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')

# Parse some parameters into vectors ----
sequencing_index_cols= unlist(strsplit(args$sequencing_index_cols, ","))
id_cols= unlist(strsplit(args$id_cols, ","))

# Validation: Check that sequencing_index_cols present in the sample meta ----
if(!validate_columns_exist(sequencing_index_cols, sample_meta)) {
  stop('One or more sequencing_index_cols is NOT present in the sample meta.')
}

# Validation: Check that id_cols are present in the sample meta ----
if(!validate_columns_exist(id_cols, sample_meta)) {
  stop('One or more id_cols is NOT present in the sample meta.')
}

# Run collate_fastq_reads on chunks of raw_counts_uncollapsed.csv ----
# raw_counts_uncollapsed can be too large to read into memory,
# so collate_fastq_reads is performed on chunks of the large file.
chunked_results= process_in_chunks(large_file_path= args$raw_counts_uncollapsed, 
                                   chunk_size= 10^6, 
                                   action= collate_fastq_reads,
                                   # Parameters for collate_fastq_reads
                                   sample_meta= sample_meta, 
                                   sequencing_index_cols= sequencing_index_cols,
                                   id_cols= id_cols,
                                   known_barcodes= unique(c(cell_line_meta$Sequence, CB_meta$Sequence)),
                                   reverse_index2= args$reverse_index2,
                                   barcode_col= args$barcode_col,
                                   low_abundance_threshold= as.numeric(args$low_abundance_threshold))

# From each chunk, extract prism_barcode_counts and bind the rows together into one dataframe.
prism_barcode_counts= data.table::rbindlist(lapply(chunked_results, function(x) x$prism_barcode_counts))
# Use data.table to group_by id_cols and barcode_col to sum up reads across all chunks.
prism_barcode_counts= prism_barcode_counts[, .(n= sum(n)), by= c(id_cols, args$barcode_col)]

# From each chunk, extract unknown_barcode_counts and bind the rows together into one dataframe.
unknown_barcode_counts= data.table::rbindlist(lapply(chunked_results, function(x) x$unknown_barcode_counts))
# Use data.table to group_by id_cols and barcode_col to sum up reads across all chunks.
unknown_barcode_counts= unknown_barcode_counts[, .(n= sum(n)), by= c(id_cols, args$barcode_col)]

# Validation: Basic file size check ----
if(nrow(prism_barcode_counts) == 0) {
  stop('ERROR: Empty file generated. No rows in prism_barcode_counts output.')
} 

# Write out files ----
out_file= paste(args$out, 'prism_barcode_counts.csv', sep='/')
print(paste("Writing prism_barcode_counts.csv to ", out_file))
write.csv(prism_barcode_counts, out_file, row.names= FALSE, quote= FALSE)

out_file= paste(args$out, 'unknown_barcode_counts.csv', sep='/')
print(paste("Writing unknown_barcode_counts.csv to ", out_file))
write.csv(unknown_barcode_counts, out_file, row.names= FALSE, quote= FALSE)
