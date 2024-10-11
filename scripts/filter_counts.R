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
source("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument('--prism_barcode_counts', default= 'prism_barcode_counts.csv', help= 'Path to prism_barcode_counts.csv')
parser$add_argument('--sample_meta', default= 'sample_meta.csv', help= 'Path to sample_meta.csv')
parser$add_argument('--cell_set_meta', default= 'cell_set_meta.csv', help= 'Path to cell_set_meta.csv')
parser$add_argument('--cell_line_meta', default= 'cell_line_meta.csv', help= 'Path to cell_line_meta.csv')
parser$add_argument("--assay_pool_meta", default="assay_pool_meta.txt", help = "Assay pool metadata")
parser$add_argument('--CB_meta', default= 'CB_meta.csv', help= 'Path to CB_meta.csv')
parser$add_argument('--id_cols', default= 'pcr_plate,pcr_well', 
                    help= 'List of sample_meta column names used to identify every PCR well')
parser$add_argument("--rm_data", type="logical", help = "Remove bad experimental data")
parser$add_argument("--pool_id", type="logical", help = "Pull pool IDs from CellDB.")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

# Read in all input files ----
prism_barcode_counts= data.table::fread(args$prism_barcode_counts, header= TRUE, sep= ',')
sample_meta= data.table::fread(args$sample_meta, header= TRUE, sep= ',')
cell_set_meta= data.table::fread(args$cell_set_meta, header= TRUE, sep= ',')
cell_line_meta= data.table::fread(args$cell_line_meta, header= TRUE, sep= ',')
CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')

# Convert input strings into vectors ----
id_cols= unlist(strsplit(args$id_cols, ","))

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

# Run filter_raw_reads -----
print('Calling filter_raw_reads ...')
module_outputs= filter_raw_reads(prism_barcode_counts= prism_barcode_counts, 
                                 sample_meta= sample_meta,
                                 cell_set_meta= cell_set_meta,
                                 cell_line_meta= cell_line_meta,
                                 CB_meta= CB_meta,
                                 id_cols= id_cols)

# Pulling pool_id when db_flag and pool_id flags are passed ----
if (args$pool_id) {
  assay_pool_meta = read.delim(args$assay_pool_meta)
  unique_cell_sets <- unique(sample_meta$cell_set[sample_meta$cell_set != ""])
  assay_pool_meta <- assay_pool_meta[assay_pool_meta$davepool_id %in% unique_cell_sets,] %>% 
    select(pool_id, ccle_name, davepool_id, depmap_id)
  
  module_outputs$filtered_counts = module_outputs$filtered_counts %>% 
    merge(assay_pool_meta, by.x=c("ccle_name", "cell_set", "depmap_id"), 
          by.y=c("ccle_name", "davepool_id", "depmap_id"), all.x=T) 
  
  module_outputs$annotated_counts = module_outputs$annotated_counts %>% 
    merge(assay_pool_meta, by.x=c("ccle_name", "cell_set", "depmap_id"), 
          by.y=c("ccle_name", "davepool_id", "depmap_id"), all.x=T)
}

# Validation: Basic file size check ----
if(sum(module_outputs$filtered_counts$n) == 0) {
  stop('All entries in filtered counts are missing!')
}

filtered_counts= module_outputs$filtered_counts
cl_entries= filtered_counts %>% dplyr::filter(!is.na(depmap_id))
if(sum(cl_entries$n) == 0) {
  stop('All cell line counts are zero!')
}
# Remove data ----
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

# Write out files ----
annot_out_file= paste0(args$out, '/annotated_counts.csv')
print(paste('Writing annotated counts to: ', annot_out_file))
module_outputs$annotated_counts %>% write.csv(annot_out_file, row.names= FALSE)

filtrc_out_file = paste(args$out, 'filtered_counts.csv', sep='/')
print(paste("Writing filtered counts csv to: ", filtrc_out_file))
write.csv(filtered_counts, filtrc_out_file, row.names=F, quote=F)
