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
parser$add_argument("--cell_line_cols", default = "depmap_id,pool_id")
parser$add_argument("--id_cols", default = "pcr_plate,pcr_well")
parser$add_argument("--count_threshold", default = 40)
parser$add_argument("--pseudocount", default = 20)
parser$add_argument("--filter_qc_flags", default = "true",
                    help = "Filter out wells with QC flags. Default is TRUE")

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

# Check if the output directory exists, if not create it
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
filter_qc_flags <- as.logical(toupper(args$filter_qc_flags))

# BY WELL (PCR_PLATE, PCR_WELL) ----------
id_cols_table <- generate_id_cols_table(
    normalized_counts = normalized_counts, annotated_counts = annotated_counts, unknown_counts = unknown_counts,
    cell_set_meta = cell_set_meta, id_cols_list = id_cols_list, cell_line_cols = cell_line_cols_list,
    count_threshold = count_threshold, cb_meta = cb_meta, pseudocount = pseudocount
)

# ID_COLS QC FLAGS
# Generate QC flags for the id_cols table and filter out flagged wells
result <- id_cols_qc_flags(
    annotated_counts = annotated_counts,
    normalized_counts = normalized_counts,
    unknown_counts = unknown_counts,
    cb_meta = cb_meta
)
id_cols_qc_flags_table <- result$well_flags
filtered_normalized_counts <- result$result

# CELL LINE BY PLATE (pcr_plate,depmap_id) ---------
# Filter out control barcodes (where depmap_id is NA) from normalized and filtered counts
normalized_counts_rm_ctl <- filtered_normalized_counts %>%
    filter_control_barcodes()
filtered_counts_rm_ctl <- filtered_counts %>%
    filter_control_barcodes()
plate_cell_table <- generate_cell_plate_table(
    normalized_counts = normalized_counts_rm_ctl, filtered_counts = filtered_counts_rm_ctl,
    cell_line_cols = cell_plate_list, pseudocount = pseudocount)

# CELL LINE PLATE QC FLAGS
# Generate QC flags for the plate_cell table and filter out flagged wells
result <- pool_well_qc_flags(
    normalized_counts = filtered_normalized_counts
)
pool_well_qc_flags_table <- result$pool_well_flags
filtered_normalized_counts <- result$result


# WRITE OUT RESULTS ----
# Write plate_cell table ----------
plate_cell_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_table.csv")
print(paste0("Writing out plate_cell_qc_table to ", plate_cell_outpath))
write.csv(
    x = plate_cell_table, file = plate_cell_outpath, row.names = FALSE,
    quote = FALSE
)
check_file_exists(plate_cell_outpath)

# Write pool_well_qc_flags table ----------
pool_well_qc_flags_outpath <- paste0(args$out, "/qc_tables/pool_well_qc_flags.csv")
print(paste0("Writing out pool_well_qc_flags to ", pool_well_qc_flags_outpath))
write.csv(
    x = pool_well_qc_flags_table, file = pool_well_qc_flags_outpath, row.names = FALSE,
    quote = FALSE
)
check_file_exists(pool_well_qc_flags_outpath)

# Write id_cols table ----------
id_cols_outpath <- paste0(args$out, "/qc_tables/id_cols_qc_table.csv")
print(paste0("Writing out id_cols_qc_table to ", id_cols_outpath))
write.csv(
    x = id_cols_table, file = id_cols_outpath, row.names = FALSE, quote = FALSE
)
check_file_exists(id_cols_outpath)

# Write id_cols_qc_flags table ----------
id_cols_qc_flags_outpath <- paste0(args$out, "/qc_tables/id_cols_qc_flags.csv")
print(paste0("Writing out id_cols_qc_flags to ", id_cols_qc_flags_outpath))
write.csv(
    x = id_cols_qc_flags_table, file = id_cols_qc_flags_outpath, row.names = FALSE,
    quote = FALSE
)
check_file_exists(id_cols_qc_flags_outpath)

if (args$filter_qc_flags) {
    # Filter out wells with QC flags
    print("Filtering out wells with QC flags")
    # Write original normalized counts ----------
    normalized_counts_original_outpath <- paste0(args$out, "/normalized_counts_original.csv")
    print(paste0("Writing unfiltered normalized_counts to ", normalized_counts_original_outpath))
    write.csv(
        x = normalized_counts, file = normalized_counts_original_outpath, row.names = FALSE,
        quote = FALSE)
    check_file_exists(normalized_counts_original_outpath)

    # Write filtered normalized counts ----------
    filtered_normalized_counts_outpath <- paste0(args$out, "/normalized_counts.csv")
    print(paste0("Writing filtered normalized_counts to ", filtered_normalized_counts_outpath))
    write.csv(
        x = filtered_normalized_counts, file = filtered_normalized_counts_outpath, row.names = FALSE,
        quote = FALSE)
    check_file_exists(filtered_normalized_counts_outpath)
    } else {
    print("Nomalized counts not filtered for qc_flags.")
    }


paste0("QC module completed.")
