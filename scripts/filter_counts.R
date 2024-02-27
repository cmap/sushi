suppressPackageStartupMessages(library(argparse))

source("./src/cellDB_metadata.R")
source("./src/filter_raw_reads.R")

suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr)) #write_delim
suppressPackageStartupMessages(library(stringr)) #str_detect
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(tidyr)) #pivot_wider
# library(prismSeqR)
suppressPackageStartupMessages(library(sets))
suppressPackageStartupMessages(library(tidyverse)) # load last - after dplyr
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(jsonlite))


## writes configuration to file
##
## takes:
##      args: args object from argparse
print_args <- function(args){
  config <- data.frame(args=names(args), values=unname(unlist(args)))
  config_path = paste(
    args$out,
    "config.txt",
    sep="/"
  )
  print(paste("Saving config.txt file in :", config_path))
  write_delim(config, config_path, delim = ": ", col_names=F)
}

# create parser object
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
parser$add_argument("--cell_line_meta", default="../metadata/cell_line_meta.csv", help = "Cell Line metadata")
parser$add_argument("--cell_set_meta", default="../metadata/cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")
parser$add_argument("--id_cols", default="cell_set,treatment,dose,dose_unit,day,bio_rep,tech_rep",
                    help = "Columns used to generate profile ids, comma-separated colnames from --sample_meta")
parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
parser$add_argument("--reverse_index2", action="store_true", default=FALSE, help = "Reverse complement of index 2 for NovaSeq and NextSeq")
parser$add_argument("--api_url", default="https://dev-api.clue.io/api/", help = "Default API URL to CellDB is DEV")
parser$add_argument("--api_key", default="", help = "Clue API key")
parser$add_argument("--db_flag", action="store_true", default=FALSE, help = "Use CellDB to locate cell set information")
parser$add_argument("--rm_data", action="store_true", default=FALSE, help = "Remove bad experimental data")
############# add a remove data flag

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

#Save arguments
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

CB_meta = read.csv(args$CB_meta)
sample_meta = read.csv(args$sample_meta)
raw_counts = read.csv(args$raw_counts)

# Using CellDB, otherwise checking static files
if (args$db_flag) {
  api_url <- args$api_url
  if (args$api_key != ""){
    api_key = args$api_key
  } else if (Sys.getenv("API_KEY") != "") {
    api_key = Sys.getenv("API_KEY")
  } else {
    stop("No API key provided via argument or environment variable API_KEY.")
  }
  
  print("Using CellDB to locate cell information.")
  cell_sets_df <- get_cell_api_info(paste(api_url,"cell_sets", sep = "/"), api_key)
  cell_pools_df <- get_cell_api_info(paste(api_url,"assay_pools", sep = "/"), api_key)
  cell_lines_df <- get_cell_api_info(paste(api_url,"cell_lines", sep = "/"), api_key)
  assay_pools_df <- get_cell_api_info(paste(api_url,"cell_set_definition_files", sep = "/"), api_key)
  
  # Renaming assay pool dataframe to act as cell_line_meta + matching case sensitivity of columns to that of static files
  cell_line_cols= c('DepMap_ID', 'CCLE_name', 'Sequence', 'LUA')
  cell_line_meta <- cell_lines_df %>%
    rename("LUA" = "lua",
           "Sequence" = "dna_sequence",
           "DepMap_ID" = "depmap_id",
           "CCLE_name" = "ccle_name") %>% dplyr::select(any_of(c(cell_line_cols)))
  
  cell_sets <- create_cell_set_meta(sample_meta, cell_sets_df, cell_pools_df)
  cell_set_meta <- cell_sets[[1]]
  failed_cell_sets <- cell_sets[[2]]
  
  # Fail if failed cell sets not empty
  if (length(failed_cell_sets) != 0) {
    print(failed_cell_sets)
    stop("The sample_meta contains the above cell sets which are not registered in CellDB:")
  }
  
  # Writing out cell_line_meta - may be necessary for downstream SUSHI?
  cell_line_out_file = paste(args$out, 'cell_line_meta.csv', sep='/')
  print(paste("writing cell_line_meta to: ", cell_line_out_file))
  write.csv(cell_line_meta, cell_line_out_file, row.names=F, quote=F)
  
  # Writing out cell_set_meta for record keeping
  cell_set_out_file = paste(args$out, 'cell_set_meta.csv', sep='/')
  print(paste("writing cell_set_meta to: ", cell_set_out_file))
  write.csv(cell_set_meta, cell_set_out_file, row.names=F, quote=F)
} else {
  print("Using static cell set information files to locate cell information.")
  cell_line_meta = read.csv(args$cell_line_meta)
  cell_set_meta = read.csv(args$cell_set_meta)
}

#split id_cols args
id_cols = unlist(strsplit(args$id_cols, ","))

if (!all(id_cols %in% colnames(sample_meta))){
  stop(paste("Colnames not found in sample_meta, check metadata or --id_cols argument:", args$id_cols))
}

sample_meta$profile_id = do.call(paste,c(sample_meta[id_cols], sep=':'))

count_threshold_arg= args$count_threshold
count_threshold = as.numeric(count_threshold_arg)

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

print("creating filtered count file")
filtered_counts = filter_raw_reads(
  raw_counts,
  sample_meta,
  cell_line_meta,
  cell_set_meta,
  CB_meta,
  id_cols=id_cols,
  count_threshold=count_threshold,
  reverse_index2=args$reverse_index2
)

# Pulling pool_id and culture when db_flag is passed
if (args$db_flag) {
  unique_cell_sets <- unique(sample_meta$cell_set[sample_meta$cell_set != ""])
  assay_pools_df <- assay_pools_df[assay_pools_df$davepool_id %in% unique_cell_sets,] %>% 
    select(pool_id, ccle_name, davepool_id, depmap_id)
  
  filtered_counts$filtered_counts = filtered_counts$filtered_counts %>% 
    merge(assay_pools_df, by.x=c("CCLE_name", "cell_set", "DepMap_ID"), by.y=c("ccle_name", "davepool_id", "depmap_id"), all.x=T) 
  
  filtered_counts$annotated_counts = filtered_counts$annotated_counts %>% 
    merge(assay_pools_df, by.x=c("CCLE_name", "cell_set", "DepMap_ID"), by.y=c("ccle_name", "davepool_id", "depmap_id"), all.x=T)
}

# Write out module outputs
qc_table = filtered_counts$qc_table
qc_out_file = paste(args$out, 'QC_table.csv', sep='/')
print(paste("writing QC_table to: ", qc_out_file))
write.csv(qc_table, qc_out_file, row.names=F, quote=F)

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
