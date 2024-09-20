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

# Validation: Check that sequencing_index_cols are from sample meta column names ----
if(!all(sequencing_index_cols %in% colnames(sample_meta))) {
  stop(paste('The following sequencing_index_cols were not found in the sample meta: ',
             sequencing_index_cols[!sequencing_index_cols %in% colnames(sample_meta)]))
}

# Validation: Check that id_cols are from sample meta column names ----
if(!all(id_cols %in% colnames(sample_meta))) {
  stop(paste('The following id_cols were not found in the sample meta: ',
             id_cols[!id_cols %in% colnames(sample_meta)]))
}

# Run collate_fastq_reads on chunks ----
# Set up loop to process chunks
header_col_names= data.table::fread(args$raw_counts_uncollapsed, header=T, sep= ',', nrow= 0) %>% colnames()
chunk_size= 10^6 # Maximum number of rows in a chunk
chunk_idx= 1 # Counter to keep track of chunks in a loop
current_chunk_size= chunk_size # Variable for loop exit condition
chunk_collector= list() # List to collect processed chunks

# For each chunk, call collate
while(current_chunk_size == chunk_size) {
  nori_chunk= data.table::fread(args$raw_counts_uncollapsed, header= F, sep= ',',
                                col.names= header_col_names,
                                nrow= chunk_size, skip= chunk_size * (chunk_idx - 1) + 1)
  
  current_chunk_size= nrow(nori_chunk) # set current chunk size to stop loop
  print(paste('Working on chunk', chunk_idx, 'with', current_chunk_size, 'rows.', sep= ' '))
  
  chunk_collector[[chunk_idx]]= collate_fastq_reads(nori_chunk, sample_meta,
                                                    sequencing_index_cols= sequencing_index_cols,
                                                    id_cols= id_cols,
                                                    reverse_index2=  args$reverse_index2,
                                                    barcode_col= args$barcode_col)
  
  chunk_idx= chunk_idx + 1
}

raw_counts= data.table::rbindlist(chunk_collector)
raw_counts= raw_counts[, .(n= sum(n)), by= c(id_cols, args$barcode_col)]

# Validation: Basic file size check ----
if(nrow(raw_counts) == 0) {
  stop('ERROR: Empty file generated. No rows in raw_counts output.')
} 

# Write out file ----
rc_out_file= paste(args$out, 'raw_counts.csv', sep='/')
print(paste("Writing raw_counts.csv to ", rc_out_file))
write.csv(raw_counts, rc_out_file, row.names= FALSE, quote= FALSE)
