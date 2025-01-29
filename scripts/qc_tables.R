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
parser$add_argument("-n", "--negcon_type", default = "ctl_vehicle")
parser$add_argument("-p", "--poscon_type", default = "trt_poscon")
parser$add_argument("--cell_line_cols", default = "depmap_id,pool_id")
parser$add_argument("--id_cols", default = "pcr_plate,pcr_well")
parser$add_argument("--count_threshold", default = 40)
parser$add_argument("--pseudocount", default = 20)

# QC Thresholds -----
# wells
parser$add_argument("--cb-tc_threshold", default = 0.5, help = "Threshold for centered cb_intercep + log2_total_count; positive and negative control wells")
parser$add_argument("--cb-mae_threshold", default = 1, help = "Threshold for mean absolute error of control barcodes; all wells")
parser$add_argument("--cb-spearman_threshold", default = 0.88, help = "Threshold for spearman correlation of control barcodes; all wells")
# pool/well combinations
parser$add_argument("--nc-pool-well-via_threshold", default = 2, help = "Remove pool/well combinations if more than half of their cell lines have < 50% or > 200% viability; negative control wells")
parser$add_argument("--pc-pool-well-via_threshold", default = 0.25, help = "Remove pool/well combinations if more than half of their cell lines have > 25%; positive control wells")
parser$add_argument("--well-pool-rm_threshold", default = (1/3), help = "If more than 1/3 of pools in a well are removed, remove the well; all wells")
# pool/pcr_plate combinations
parser$add_argument("--ctl-well-plate_threshold", default = (1/3), help = "If more than 1/3 of the control wells are removed for a given pool, remove all remaining wells in that PCR plate")
# cell_line/pcr_plate combinations
parser$add_argument("--nc-cell-plate-mad_threshold", default = 1, help = "If the median absolute deviation of the negative control cell lines in a pcr_plate is greater than this value, remove the cell line from that pcr_plate")
parser$add_argument("--cell-plate-er_threshold", default = 0.05, help = "If the error rate of a cell line in a pcr_plate is greater than this value, remove the cell line from that pcr_plate")
parser$add_argument("--pc-cell-plate-via_threshold", default = 0.25, help = "If the viability of the positive controls for a given cell line in a pcr_plate is greater than this value, remove the cell line from that pcr_plate")



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
plate_cell_table <- generate_cell_plate_table(
    normalized_counts = normalized_counts, filtered_counts = filtered_counts,
    cell_line_cols = cell_plate_list, pseudocount = pseudocount)

# Write to file ----------
plate_cell_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_table.csv")
print(paste0("Writing out plate_cell_qc_table to ", plate_cell_outpath))
write.csv(
    x = plate_cell_table, file = plate_cell_outpath, row.names = FALSE,
    quote = FALSE
)
check_file_exists(plate_cell_outpath)


# BY WELL (PCR_PLATE, PCR_WELL) ----------
id_cols_table <- generate_id_cols_table(
    normalized_counts = normalized_counts, annotated_counts = annotated_counts,
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
