options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
library(jsonlite)
source ("./src/qc_functions.R")
source ("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
parser$add_argument("--cell_set_and_pool_meta", default="cell_set_and_pool_meta.csv", help= "Cell line metadata")
parser$add_argument("--normalized_counts", default= "normalized_counts.csv", help= "normalized counts file")
parser$add_argument("--annotated_counts", default= "annotated_counts.csv", help= "annotated counts file")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")
parser$add_argument("-n", "--negcon_type", default= "ctl_vehicle")
parser$add_argument("-p", "--poscon_type", default= "trt_poscon")

args <- parser$parse_args()

# Read in metadata files as data.table objects ----
cell_set_meta= data.table::fread(args$cell_set_and_pool_meta, header= TRUE, sep= ',')
normalized_counts= data.table::fread(args$normalized_counts, header= TRUE, sep= ',')
annotated_counts= data.table::fread(args$annotated_counts, header= TRUE, sep= ',')

# CELL LINE BY PLATE (pcr_plate,depmap_id) ----------

cell_line_plate_grouping <- c("depmap_id", "pcr_plate") # Define columns to group by

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
  negcon = "ctl_vehicle",
  poscon = "trt_poscon"
)

# Merge and compute poscon LFC
plate_cell_table <- medians_and_mad %>%
  left_join(error_rates, by = group_cols) %>%
  compute_control_lfc(negcon = args$negcon_type, poscon = args$poscon_type)

# Write to file ----------
plate_cell_outpath <- paste0(args$out, "/qc_tables/plate_cell_qc_table.csv")
print(paste0("Writing out plate_cell_qc_table to ", plate_cell_outpath))
write.csv(x = plate_cell_table, file = plate_cell_outpath, row.names = FALSE, quote = FALSE)
check_file_exists(plate_cell_outpath)

# BY ID_COLS (PCR_PLATE, PCR_WELL) ----------

id_cols_grouping = c("pcr_plate", "pcr_well") # Define columns to group by
id_cols_table <- create_id_cols_table(annotated_counts = annotated_counts, group_cols = id_cols_grouping,
                                      cell_set_meta = cell_set_meta, metric = 'n')

# Write to file ----------
id_cols_outpath <- paste0(args$out, "/qc_tables/id_cols_qc_table.csv")
print(paste0("Writing out id_cols_qc_table to ", id_cols_outpath))
write.csv(x = id_cols_table, file = id_cols_outpath, row.names = FALSE, quote = FALSE)
check_file_exists(id_cols_outpath)