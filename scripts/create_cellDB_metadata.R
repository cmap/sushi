options(cli.unicode = FALSE)
suppressPackageStartupMessages(library(argparse))
source("./src/cellDB_metadata.R")
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
cell_sets_df <- get_cell_api_info(paste(api_url,"cell_sets", sep = "/"), api_key)
cell_pools_df <- get_cell_api_info(paste(api_url,"assay_pools", sep = "/"), api_key)
cell_lines_df <- get_cell_api_info(paste(api_url,"cell_lines", sep = "/"), api_key)
assay_pools_df <- get_cell_api_info(paste(api_url,"cell_set_definition_files", sep = "/"), api_key)
assay_pools_meta <- select(assay_pools_df, -cell_set_desc)
if (cb_ladder != "cb_custom.csv"){
  control_bc_df <- get_cell_api_info(paste(api_url,"v_control_barcodes", sep = "/"), api_key, filter = list(where = list(set = cb_ladder), fields = c("sequence", "name", "log_dose")))
} else {
  file_path <- file.path(args$out, cb_ladder)
  control_bc_df <- read.csv(file_path)
}

# Renaming assay pool dataframe to act as cell_line_meta + matching case sensitivity of columns to that of static files
cell_line_cols= c('DepMap_ID', 'ccle_name', 'Sequence', 'LUA')
cell_line_meta <- cell_lines_df %>%
  rename("LUA" = "lua",
         "Sequence" = "dna_sequence") %>% dplyr::select(any_of(c(cell_line_cols)))

# Checking if the selected cb_ladder returned any data + adjusting case sensitivity of the headers to match the original CB_meta static file
if (nrow(control_bc_df) > 0 & cb_ladder != "cb_custom.csv") {
  CB_meta <- control_bc_df %>% rename("Sequence" = "sequence", "Name" = "name")
} else {
  print(paste("Since the cb_ladder selected was '", cb_ladder, "', no renaming is necessary."))
  CB_meta <- control_bc_df
}

cell_sets <- create_cell_set_meta(sample_meta, cell_sets_df, cell_pools_df, cell_line_meta)
cell_set_meta <- cell_sets[[1]]
failed_cell_sets <- cell_sets[[2]]

# Fail if failed cell sets not empty
if (length(failed_cell_sets) != 0) {
  print(failed_cell_sets)
  stop("The sample_meta contains the above cell sets which are not registered in CellDB:")
}

# Writing out cell_line_meta
cell_line_out_file = paste(args$out, 'cell_line_meta.csv', sep='/')
print(paste("writing cell_line_meta to: ", cell_line_out_file))
write.csv(cell_line_meta, cell_line_out_file, row.names=F, quote=F)

# Writing out cell_set_meta
cell_set_out_file = paste(args$out, 'cell_set_meta.csv', sep='/')
print(paste("writing cell_set_meta to: ", cell_set_out_file))
write.csv(cell_set_meta, cell_set_out_file, row.names=F, quote=F)

# Writing out assay_pool_df
assay_pool_out_file = paste(args$out, 'assay_pool_meta.txt', sep='/')
print(paste("writing assay_pools to: ", assay_pool_out_file))
write.table(assay_pools_meta, assay_pool_out_file, row.names=F, quote=F, sep="\t")

# Writing out control_barcode_df if cb_ladder returned data in control_bc_df
if (nrow(control_bc_df) > 0) {
  control_barcode_out_file = paste(args$out, 'CB_meta.csv', sep='/')
  print(paste("writing CB_meta to: ", control_barcode_out_file))
  write.csv(CB_meta, control_barcode_out_file, row.names=F, quote=F)
} 
