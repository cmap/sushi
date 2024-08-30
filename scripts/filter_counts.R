options(cli.unicode = FALSE)
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

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--raw_counts", default="raw_counts.csv", help = "path to file containing raw counts")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("--cell_line_meta", default="cell_line_meta.csv", help= "Cell Line metadata")
parser$add_argument("--cell_set_meta", default="cell_set_meta.csv", help= "Cell set metadata")
parser$add_argument("--assay_pool_meta", default="assay_pool_meta.txt", help = "Assay pool metadata")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")
parser$add_argument("--id_cols", default= "pcr_plate,pcr_well", 
                    help = "Sequencing columns in the sample meta")
parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
parser$add_argument("--rm_data", type="logical", help = "Remove bad experimental data")
parser$add_argument("--pool_id", type="logical", help = "Pull pool IDs from CellDB.")
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
id_cols= unlist(strsplit(args$id_cols, ","))
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
filtered_counts = filter_raw_reads(raw_counts= raw_counts, sample_meta= sample_meta,
                                   cell_line_meta= cell_line_meta,
                                   cell_set_meta= cell_set_meta,
                                   CB_meta= CB_meta,
                                   id_cols= id_cols,
                                   count_threshold= as.numeric(args$count_threshold))

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
unmapped_reads= filtered_counts$unmapped_reads
unmapped_out = paste(args$out, 'unmapped_reads.csv', sep='/')
print(paste("writing unmapped reads to: ", unmapped_out))
write.csv(unmapped_reads, unmapped_out, row.names=F)

annotated_counts = filtered_counts$annotated_counts
annot_out_file = paste(args$out, 'annotated_counts.csv', sep='/')
print(paste("writing annotated counts to: ", annot_out_file))
write.csv(annotated_counts, annot_out_file, row.names=F)

filtered_counts = filtered_counts$filtered_counts

print(paste("rm_data:", args$rm_data))
# Remove data if needed
if(args$rm_data == TRUE){
  print('rm_data is TRUE, removing supplied data.')
  data_to_remove <- read.csv(paste(args$out, 'data_to_remove.csv', sep='/'))
  print('Data to remove:')
  print(head(data_to_remove))
  filt_rm <- remove_data(filtered_counts, data_to_remove)
  
  # keep the full filtered counts with the data that needs to be removed
  filtered_counts_original <- filtered_counts
  write.csv(filtered_counts_original, paste(args$out, 'filtered_counts_original.csv', sep='/'), row.names=F, quote=F)
  # re-point to what filtered_counts should be
  filtered_counts <- filt_rm
  rows_removed = nrow(filtered_counts_original) - nrow(filtered_counts)
  paste("Number of rows removed: ", rows_removed)
}

filtrc_out_file = paste(args$out, 'filtered_counts.csv', sep='/')
print(paste("writing filtered counts csv to: ", filtrc_out_file))
write.csv(filtered_counts, filtrc_out_file, row.names=F, quote=F)

