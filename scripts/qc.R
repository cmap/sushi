options(cli.unicode = FALSE)
suppressPackageStartupMessages({
  library(argparse)
  library(magrittr)
  library(tidyverse)
  library(data.table)
  library(jsonlite)
})
source ("./src/qc_functions.R")
source ("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
parser$add_argument("--cell_set_and_pool_meta", default="cell_set_and_pool_meta.csv", help= "Cell line metadata")
parser$add_argument("--normalized_counts", default= "normalized_counts.csv", help= "normalized counts file")
parser$add_argument("--annotated_counts", default= "annotated_counts.csv", help= "annotated counts file")
parser$add_argument("--filtered_counts", default= "filtered_counts.csv", help= "filtered counts file")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")
parser$add_argument("-n", "--negcon_type", default= "ctl_vehicle")
parser$add_argument("-p", "--poscon_type", default= "trt_poscon")
parser$add_argument("--cell_line_cols", default= "depmap_id,pool_id")
parser$add_argument("--id_cols", default= "pcr_plate,pcr_well")

args <- parser$parse_args()

# Read in metadata files as data.table objects ----
paste0("Reading in ", args$cell_set_and_pool_meta, ".....")
cell_set_meta <- data.table::fread(args$cell_set_and_pool_meta, header= TRUE, sep= ',')
paste0("Reading in ", args$normalized_counts, ".....")
normalized_counts <- data.table::fread(args$normalized_counts, header= TRUE, sep= ',')
paste0("Reading in ", args$annotated_counts, ".....")
annotated_counts <- data.table::fread(args$annotated_counts, header= TRUE, sep= ',')
paste0("Reading in ", args$filtered_counts, ".....")
filtered_counts <- data.table::fread(args$filtered_counts, header= TRUE, sep= ',')

# Create qc_table output directory if it doesn't exist ----
paste0("Creating output directory ", args$out, "/qc_tables.....")
if (!dir.exists(paste0(args$out, "/qc_tables"))) {
  dir.create(paste0(args$out, "/qc_tables"))
}

# CELL LINE BY PLATE (pcr_plate,depmap_id) ---------

cell_line_list <- strsplit(cell_line_cols, ",")[[1]]
cell_line_plate_grouping <- c(cell_line_list,"pcr_plate") # Define columns to group by
paste0("Computing QC metrics grouping by ", paste0(cell_line_plate_grouping, collapse = ","), ".....")

# Compute control medians and MAD
medians_and_mad <- compute_ctl_medians_and_mad(
  df = normalized_counts,
  group_cols = cell_line_plate_grouping,
  negcon = "ctl_vehicle",
  poscon = "trt_poscon"
)

# Compute error rate
error_rates <- compute_error_rate(
  df = normalized_counts,
  metric = "log2_normalized_n",
  group_cols = cell_line_plate_grouping,
  negcon = args$negcon_type,
  poscon = args$poscon_type
)

# Compute poscon LFC
poscon_lfc <- compute_control_lfc(
  df = medians_and_mad,
  negcon = args$negcon_type,
  poscon = args$poscon_type
)

# Compute cell line fractions per plate
cell_line_fractions <- compute_cl_fractions(
  df = filtered_counts,
  grouping_cols = cell_line_plate_grouping
)

# Merge all tables together
paste0("Merging ", paste0(cell_line_plate_grouping, collapse = ","), " QC tables together.....")
plate_cell_table <- medians_and_mad %>%
  dplyr::left_join(error_rates, by = cell_line_plate_grouping) %>%
  dplyr::left_join(poscon_lfc, by = cell_line_plate_grouping) %>%
  dplyr::left_join(cell_line_fractions, by = cell_line_plate_grouping)

# Write to file ----------
plate_cell_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_table.csv")
print(paste0("Writing out plate_cell_qc_table to ", plate_cell_outpath))
write.csv(x = plate_cell_table, file = plate_cell_outpath, row.names = FALSE, quote = FALSE)
check_file_exists(plate_cell_outpath)

# BY WELL (PCR_PLATE, PCR_WELL) ----------

id_cols_list <- strsplit(id_cols, ",")[[1]]
paste0("Computing QC metrics grouping by ", paste0(id_cols_list, collapse = ","), ".....")

read_stats <- compute_read_stats(annotated_counts = annotated_counts, group_cols = id_cols_list,
                                 cell_set_meta= cell_set_meta, metric = "n")

skew <- compute_skew(annotated_counts, group_cols = id_cols_list, metric = "n")

id_cols_table <- read_stats %>%
  dplyr::left_join(skew, by = id_cols_list)

# Write to file ----------
paste0("Merging ", paste0(id_cols_list, collapse = ","), " QC tables together.....")
id_cols_outpath <- paste0(args$out, "/qc_tables/id_cols_qc_table.csv")
print(paste0("Writing out id_cols_qc_table to ", id_cols_outpath))
write.csv(x = id_cols_table, file = id_cols_outpath, row.names = FALSE, quote = FALSE)
check_file_exists(id_cols_outpath)

paste0("QC module completed.")
