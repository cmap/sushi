options(cli.unicode = FALSE)
suppressPackageStartupMessages(library(argparse))
source("create_cell_meta/cellDB_metadata.R")
source("src/kitchen_utensils.R")
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr)) #write_delim
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(sets))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("--cb_ladder", default="", help = "Control barcode ladder")
parser$add_argument("--api_url", default="https://api.clue.io/api/", help = "Default API URL to CellDB is DEV")
parser$add_argument("--api_key", default="", help = "Clue API key")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

#Save arguments
if (args$out == ""){
  args$out = args$wkdir
}

# Fail if control barcode ladder not provided 
if (args$cb_ladder != ""){
  cb_ladder <- args$cb_ladder
  print(paste("Selected cb_ladder:", cb_ladder))
} else {
  stop("Control barcode ladder not specified in arguments.")
}

sample_meta = read.csv(args$sample_meta)
api_url <- args$api_url

if (args$api_key != ""){
  api_key = args$api_key
} else if (Sys.getenv("API_KEY") != "") {
  api_key = Sys.getenv("API_KEY")
} else {
  stop("No API key provided via argument or environment variable API_KEY.")
}

print("Using CellDB to locate cell information.")
print(api_url)
cell_sets_df <- get_cell_api_info(paste(api_url,"e_cell_sets", sep = "/"), api_key)
cell_pools_df <- get_cell_api_info(paste(api_url,"e_assay_pools", sep = "/"), api_key)
cell_lines_df <- get_cell_api_info(paste(api_url,"e_cell_lines", sep = "/"), api_key)
assay_pools_df <- get_cell_api_info(paste(api_url,"e_cell_set_definition_files", sep = "/"), api_key)
assay_pools_meta <- select(assay_pools_df, -cell_set_desc)

# If a custom cb_ladder is not provided, pull from CellDB
if (!str_detect(cb_ladder, ".csv")) {
  control_bc_df <- get_cell_api_info(
    paste(api_url, "v_control_barcodes", sep = "/"),
    api_key,
    filter = list(
      where = list(set = cb_ladder)
    )
  )

  # Rename columns and adjust case sensitivity to match the original CB_meta static file
  # Select only the columns that are needed
  control_bc_df <- control_bc_df %>%
    rename(
      cb_name = name,
      forward_read_barcode = sequence,
      cb_log10_dose = log_dose,
      cb_ladder = set
    ) %>%
    mutate(cb_ladder = tolower(cb_ladder)) %>%
    select(forward_read_barcode, cb_name, cb_log10_dose, cb_ladder)

} else {
  file_path <- file.path(args$out, cb_ladder)
  control_bc_df <- read.csv(file_path)
}


# Renaming assay pool dataframe to act as cell_line_meta + matching case sensitivity of columns to that of static files
cell_line_cols= c('depmap_id', 'forward_read_barcode', 'lua')
cell_line_meta <- cell_lines_df %>%
  rename("forward_read_barcode" = "dna_sequence") %>% dplyr::select(any_of(c(cell_line_cols))) %>%
  dplyr::distinct()

cell_sets <- create_cell_set_meta(sample_meta, cell_sets_df, cell_pools_df, cell_line_meta)
cell_set_meta <- cell_sets[[1]]
failed_cell_sets <- cell_sets[[2]]

# Fail if failed cell sets not empty
if (length(failed_cell_sets) != 0) {
  print(failed_cell_sets)
  stop("The sample_meta contains the above cell sets which are not registered in CellDB:")
}

# Pivoting cell_set_meta into long form 
cell_set_meta_long <- cell_set_meta %>%
  separate_rows(members, sep = ";")

# Combining cell_set_meta_long + assay_pools_meta
# Only pull pool_id if all "cell_set" values passed in sample_meta exist in "davepool_id" of assay_pools_meta
if(all(cell_set_meta_long$cell_set %in% assay_pools_meta$davepool_id)) {
  print("Merging cell set metadata with assay pool metadata to pull pool_id.")
  cell_set_assay_pool_meta <- cell_set_meta_long %>%
    inner_join(assay_pools_meta, by = c("cell_set" = "davepool_id", "members" = "depmap_id")) %>%
    select(cell_set, pool_id, barcode_id, depmap_id = members) %>%
    rename("lua" = "barcode_id") %>%
    dplyr::distinct()
} else {
  print("One or more cell sets passed in sample_meta have not been registered in CellDB. Unable to pull pool_id in cell_set_meta.")
  cell_set_assay_pool_meta <- cell_set_meta_long %>%
    select(cell_set, depmap_id = members) %>%
    dplyr::distinct()
}

CB_meta <- control_bc_df

# Join CB_meta with cell_line_meta if the depmap_id and lua columns are not present in CB_meta
if (!all(c("depmap_id", "lua") %in% colnames(CB_meta))) {
  print("Adding depmap_id and lua columns to CB_meta.")
  CB_meta <- CB_meta %>%
    dplyr::left_join(cell_line_meta, by = "forward_read_barcode")
}

# Writing out cell_line_meta
cell_line_out_file = paste(args$out, 'cell_line_meta.csv', sep='/')
print(paste("writing cell_line_meta to: ", cell_line_out_file))
write.csv(cell_line_meta, cell_line_out_file, row.names=F, quote=F)

# Writing out combined cell_set_meta and assay_pool_meta
cell_set_assay_pool_out_file = paste(args$out, 'cell_set_and_pool_meta.csv', sep='/')
print(paste("writing combined cell_set_meta and assay_pool_meta file to: ", cell_set_assay_pool_out_file))
write.csv(cell_set_assay_pool_meta, cell_set_assay_pool_out_file, row.names=F, quote=F)

# Writing out control_barcode_df if cb_ladder returned data in control_bc_df
if (nrow(control_bc_df) > 0) {
  control_barcode_out_file = paste(args$out, 'CB_meta.csv', sep='/')
  print(paste("writing CB_meta to: ", control_barcode_out_file))
  write.csv(CB_meta, control_barcode_out_file, row.names=F, quote=F)
}

# Ensure that cell_line_meta and cell_set_assay_pool_meta were successfully generated
check_file_exists(cell_line_out_file)
check_file_exists(cell_set_assay_pool_out_file)
