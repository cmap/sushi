suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr)) #write_delim
suppressPackageStartupMessages(library(stringr)) #str_detect
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(tidyr)) #pivot_wider
suppressPackageStartupMessages(library(sets))
suppressPackageStartupMessages(library(tidyverse)) # load last - after dplyr
source("./src/filter_raw_reads.R")


## writes configuration to file
##
## takes:
##      args: args object from argparse
print_args <- function(args){
  config <- data.frame(args=names(args), values=unname(unlist(args)))
  config_path = paste(args$out, "config.txt", sep="/")
  print(paste("Saving config.txt file in :", config_path))
  write_delim(config, config_path, delim = ": ", col_names=F)
}

# Arguement parser ----
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--raw_counts", default="raw_counts.csv", help = "path to file containing raw counts")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("--cell_line_meta", default="cell_line_meta.csv", help = "Cell Line metadata")
parser$add_argument("--cell_set_meta", default="cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("--assay_pool_meta", default="assay_pool_meta.txt", help = "Assay pool metadata")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")
parser$add_argument("--sequencing_index_cols", default= "index_1,index_2", 
                    help = "Sequencing columns in the sample meta")
parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
parser$add_argument("--reverse_index2", default=FALSE,
                    help = "Reverse complement of index 2 for NovaSeq and NextSeq")
parser$add_argument("--rm_data", default=FALSE, help = "Remove bad experimental data")
parser$add_argument("--pool_id", default=FALSE, help = "Pull pool IDs from CellDB.")
parser$add_argument("--control_type", default="negcon", 
                    help = "negative control wells in trt_type column in sample metadata")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

# Read in files and set up parameters ----
cell_set_meta= data.table::fread(args$cell_set_meta, header= T, sep= ',', data.table= F)
cell_line_meta= data.table::fread(args$cell_line_meta, header= T, sep= ',', data.table= F)
CB_meta= data.table::fread(args$CB_meta, header= T, sep= ',', data.table= F)
sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',', data.table= F)
raw_counts= data.table::fread(args$raw_counts, header= T, sep= ',', data.table= F)

# Convert strings to vectors ----
# Also check that column names are present in the sample meta.
sequencing_index_cols= unlist(strsplit(args$sequencing_index_cols, ","))
if (!all(sequencing_index_cols %in% colnames(sample_meta))){
  stop(paste("All seq columns not found in sample_meta, check metadata or --sequencing_index_cols argument:",
             args$sequencing_index_cols))
}

count_threshold = as.numeric(args$count_threshold)

# make sure LUA codes in cell line meta are unique
cell_line_meta %<>% 
  dplyr::group_by(LUA) %>% 
  dplyr::mutate(LUA.duplicity = n()) %>% 
  dplyr::ungroup()

print(paste0("LUAs that are duplicated ", 
             dplyr::filter(cell_line_meta, LUA.duplicity > 1)$LUA %>% 
               unique() %>% sort() %>% paste0(collapse = ", "))) # print LUA duplicates

cell_line_meta %<>% 
  dplyr::filter(!duplicated(cell_line_meta$LUA, fromLast = TRUE)) %>%
  dplyr::select(-LUA.duplicity)

# Remove flowcell_name and lane columns from sample_meta because
# there is a profile_id duplicate when there are more than 1 seq runs
#sample_meta %<>% select(-flowcell_name, -flowcell_lane) %>%
 # distinct() # This needs to be removed for sequencing_index_cols to work! - YL

# Run filter_raw_reads -----
print("creating filtered count file")
filtered_counts = filter_raw_reads(raw_counts,
                                   sample_meta,
                                   cell_line_meta,
                                   cell_set_meta,
                                   CB_meta,
                                   sequencing_index_cols= sequencing_index_cols,
                                   count_threshold= as.numeric(args$count_threshold),
                                   reverse_index2= args$reverse_index2)

# Pulling pool_id when db_flag and pool_id flags are passed
if (args$pool_id) {
  assay_pool_meta = read.delim(args$assay_pool_meta)
  unique_cell_sets <- unique(sample_meta$cell_set[sample_meta$cell_set != ""])
  assay_pool_meta <- assay_pool_meta[assay_pool_meta$davepool_id %in% unique_cell_sets,] %>% 
    select(pool_id, ccle_name, davepool_id, depmap_id)
  
  filtered_counts$filtered_counts = filtered_counts$filtered_counts %>% 
    merge(assay_pool_meta, by.x=c("CCLE_name", "cell_set", "DepMap_ID"), 
          by.y=c("ccle_name", "davepool_id", "depmap_id"), all.x=T) 
  
  filtered_counts$annotated_counts = filtered_counts$annotated_counts %>% 
    merge(assay_pool_meta, by.x=c("CCLE_name", "cell_set", "DepMap_ID"), 
          by.y=c("ccle_name", "davepool_id", "depmap_id"), all.x=T)
}

# Validation: Basic file size check ----
if(sum(filtered_counts$filtered_counts$n) == 0) {
  stop('All entries in filtered counts are missing!')
}

cl_entries= filtered_counts$filtered_counts %>% dplyr::filter(!is.na(CCLE_name))
if(sum(cl_entries$n) == 0) {
  stop('All cell line counts are zero!')
}

# Write out module outputs ----
qc_table = filtered_counts$qc_table
qc_out_file = paste(args$out, 'QC_table.csv', sep='/')
print(paste("writing QC_table to: ", qc_out_file))
write.csv(qc_table, qc_out_file, row.names=F, quote=F)

unmapped_reads= filtered_counts$unmapped_reads
unmapped_out = paste(args$out, 'unmapped_reads.csv', sep='/')
print(paste("writing unmapped reads to: ", unmapped_out))
write.csv(unmapped_reads, unmapped_out, row.names=F)

annotated_counts = filtered_counts$annotated_counts
annot_out_file = paste(args$out, 'annotated_counts.csv', sep='/')
print(paste("writing annotated counts to: ", annot_out_file))
write.csv(annotated_counts, annot_out_file, row.names=F)

filtered_counts = filtered_counts$filtered_counts

# Remove data if needed
if(args$rm_data == T){
  data_to_remove <- read.csv(paste(args$out, 'data_to_remove.csv', sep='/'))
  filt_rm <- remove_data(filtered_counts, data_to_remove)
  
  # keep the full filtered counts with the data that needs to be removed
  filtered_counts_original <- filtered_counts
  write.csv(filtered_counts_original, paste(args$out, 'filtered_counts_original.csv', sep='/'), row.names=F, quote=F)
  # re-point to what filtered_counts should be
  filtered_counts <- filt_rm
}

filtrc_out_file = paste(args$out, 'filtered_counts.csv', sep='/')
print(paste("writing filtered counts csv to: ", filtrc_out_file))
write.csv(filtered_counts, filtrc_out_file, row.names=F, quote=F)

