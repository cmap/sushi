library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
library(jsonlite)
source ("/Users/jdavis/code/sushi/scripts/src/qc_functions.R")

# Argument parser ----
parser <- ArgumentParser()
parser$add_argument("--raw_counts_uncollapsed", default= "raw_counts_uncollapsed.csv",
                    help= "path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("--cell_line_meta", default="cell_line_meta.csv", help= "Cell line metadata")
parser$add_argument("--CB_meta", default= "CB_meta.csv", help= "Control Barcode metadata")
parser$add_argument("--l2fc", default= "l2fc.csv", help= "l2fc file")
parser$add_argument("--collapsed_l2fc", default= "collapsed_l2fc.csv", help= "collapsed l2fc file")
parser$add_argument("--normalized_counts", default= "normalized_counts.csv", help= "normalized counts file")
parser$add_argument("--annotated_counts", default= "annotated_counts.csv", help= "annotated counts file")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")
parser$add_argument("-c", "--config", default= "config.json", help= "Configuration file")
parser$add_argument("-sig", "--sig_cols", default= "cell_set,pert_name,pert_dose,pert_dose_unit,day,pert_type,x_project_id,pert_plate")
parser$add_argument("-id", "--id_cols", default= "pcr_plate,pcr_well")
parser$add_argument("-negcon", "--negcon_type", default= "ctl_vehicle")
parser$add_argument("-poscon", "--poscon_type", default= "trt_poscon")

parser <- parser$parse_args()

# Read in metadata files as data.table objects ----
#sample_meta= data.table::fread(args$sample_meta, header= TRUE, sep= ',')
cell_set_meta= data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/cell_set_and_pool_meta.csv", header= TRUE, sep= ',')
#CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')
#l2fc= data.table::fread(args$l2fc, header= TRUE, sep= ',')
#collapsed_l2fc= data.table::fread(args$collapsed_l2fc, header= TRUE, sep= ',')
normalized_counts= data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/normalized_counts.csv", header= TRUE, sep= ',')
annotated_counts= data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/annotated_counts.csv", header= TRUE, sep= ',')

# CELL LINE BY PLATE (pcr_plate,depmap_id) ----------
group_cols <- c("depmap_id", "pcr_plate") # Define columns to group by

# Compute control medians and MAD
medians_and_mad <- compute_ctl_medians_and_mad(
  df = normalized_counts,
  group_cols = group_cols,
  negcon = "ctl_vehicle",
  poscon = "trt_poscon"
)

# Compute error rate
error_rates <- compute_error_rate(
  df = normalized_counts,
  metric = "log2_normalized_n",
  group_cols = group_cols,
  negcon = "ctl_vehicle",
  poscon = "trt_poscon"
)

# Merge and compute poscon LFC
plate_cell_table <- medians_and_mad %>%
  left_join(error_rates, by = group_cols) %>%
  compute_control_lfc(negcon = "ctl_vehicle", poscon = "trt_poscon")

# BY ID_COLS (PCR_PLATE, PCR_WELL) ----------

id_cols_table <- create_id_cols_table(annotated_counts = annotated_counts, group_cols = c("pcr_plate", "pcr_well"),
                                      cell_set_meta = cell_set_meta, metric = 'n')