options(cli.unicode = FALSE)
suppressPackageStartupMessages(
    {
        library(argparse)
        library(magrittr)
        library(tidyverse)
        library(data.table)
        library(jsonlite)
    }
)
source("./src/qc_table_functions.R")
source("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
parser$add_argument(
    "--cell_set_and_pool_meta", default = "cell_set_and_pool_meta.csv",
    help = "Cell line metadata"
)
parser$add_argument(
    "--normalized_counts", default = "normalized_counts.csv", help = "normalized counts file"
)
parser$add_argument(
    "--annotated_counts", default = "annotated_counts.csv", help = "annotated counts file"
)
parser$add_argument(
    "--filtered_counts", default = "filtered_counts.csv", help = "filtered counts file"
)
parser$add_argument(
    "-o", "--out", default = getwd(), help = "Output path. Default is working directory"
)
parser$add_argument(
    "--control_barcode_meta", default = "CB_meta.csv", help = "Control barcode metadata"
)
parser$add_argument(
    "--unknown_barcode_counts", default = "unknown_barcode_counts.csv",
    help = "Unknown barcode counts file"
)
parser$add_argument("-n", "--negcon_type", default = "ctl_vehicle")
parser$add_argument("-p", "--poscon_type", default = "trt_poscon")
parser$add_argument("--cell_line_cols", default = "depmap_id,pool_id,lua")
parser$add_argument("--id_cols", default = "pcr_plate,pcr_well")
parser$add_argument("--count_threshold", default = 40)
parser$add_argument("--pseudocount", default = 20)

args <- parser$parse_args()

# Read in metadata files as data.table objects ----
paste0("Reading in ", args$cell_set_and_pool_meta, ".....")
cell_set_meta <- data.table::fread(args$cell_set_and_pool_meta, header = TRUE, sep = ",")
paste0("Reading in ", args$normalized_counts, ".....")
normalized_counts <- data.table::fread(args$normalized_counts, header = TRUE, sep = ",")
paste0("Reading in ", args$annotated_counts, ".....")
annotated_counts <- data.table::fread(args$annotated_counts, header = TRUE, sep = ",")
paste0("Reading in ", args$filtered_counts, ".....")
filtered_counts <- data.table::fread(args$filtered_counts, header = TRUE, sep = ",")
paste0("Reading in ", args$control_barcode_meta, ".....")
cb_meta <- data.table::fread(args$control_barcode_meta, header = TRUE, sep = ",")
paste0("Reading in ", args$unknown_barcode_counts, header = TRUE, sep = ",")
unknown_counts <- data.table::fread(args$unknown_barcode_counts, header = TRUE, sep = ",")

# Create qc_table output directory if it doesn't exist ----
paste0("Creating output directory ", args$out, "/qc_tables.....")
if (!dir.exists(paste0(args$out, "/qc_tables")))
    {
    dir.create(paste0(args$out, "/qc_tables"))
}

# DEFINE COLUMNS
cell_line_cols <- args$cell_line_cols
cell_line_cols_list <- strsplit(cell_line_cols, ",")[[1]]
cell_plate_list <- c(cell_line_cols, "pcr_plate")

id_cols <- args$id_cols
id_cols_list <- strsplit(id_cols, ",")[[1]]

count_threshold <- as.numeric(args$count_threshold)

pseudocount <- as.numeric(args$pseudocount)

# CELL LINE BY PLATE (pcr_plate,depmap_id) ---------
# Filter out control barcodes (where depmap_id is NA) from normalized and filtered counts
normalized_counts_rm_ctl <- normalized_counts %>%
    filter_control_barcodes()
filtered_counts_rm_ctl <- filtered_counts %>%
    filter_control_barcodes()
plate_cell_table <- generate_cell_plate_table(
    normalized_counts = normalized_counts_rm_ctl, filtered_counts = filtered_counts_rm_ctl,
    cell_line_cols = cell_plate_list, pseudocount = pseudocount)

# Write to file for internal use ----------
plate_cell_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_table_internal.csv")
print(paste0("Writing out internal plate_cell_qc_table to ", plate_cell_outpath))
write.csv(
    x = plate_cell_table, file = plate_cell_outpath, row.names = FALSE,
    quote = FALSE
)
check_file_exists(plate_cell_outpath)


# Write to file for portal use----------
plate_cell_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_table.csv")
print(paste0("Writing out external plate_cell_qc_table to ", plate_cell_outpath))
write.csv(
    x = plate_cell_table %>% 
        dplyr::select(
            c("pool_id", "depmap_id", "lua", "pcr_plate",
              "pert_plate",
              "error_rate", "lfc_trt_poscon",  
              "median_raw_ctl_vehicle", "mad_log_normalized_ctl_vehicle",  
              "median_log_normalized_ctl_vehicle", 
              "n_replicates_ctl_vehicle", "n_replicates_trt_poscon", 
              "viability_trt_poscon", "qc_pass", "qc_pass_pert_plate")),
    file = plate_cell_outpath, row.names = FALSE,
    quote = FALSE
)
check_file_exists(plate_cell_outpath)


# BY WELL (PCR_PLATE, PCR_WELL) ----------
id_cols_table <- generate_id_cols_table(
    normalized_counts = normalized_counts, annotated_counts = annotated_counts, unknown_counts = unknown_counts,
    cell_set_meta = cell_set_meta, id_cols_list = id_cols_list, cell_line_cols = cell_line_cols_list,
    count_threshold = count_threshold, cb_meta = cb_meta, pseudocount = pseudocount
)

# Write to file ----------
paste0(
    "Merging ", paste0(id_cols_list, collapse = ","),
    " QC tables together....."
)
id_cols_outpath <- paste0(args$out, "/qc_tables/id_cols_qc_table.csv")
print(paste0("Writing out id_cols_qc_table to ", id_cols_outpath))
write.csv(
    x = id_cols_table, file = id_cols_outpath, row.names = FALSE, quote = FALSE
)
check_file_exists(id_cols_outpath)

paste0("QC module completed.")
