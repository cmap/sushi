library(argparse)
library(tidyverse)

source("utils/kitchen_utensils.R")
source("bias_correction/bias_correction_functions.R")

# Shell script argument parser ----
parser = ArgumentParser()

parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE, help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose", help = "Print little output")
parser$add_argument("--l2fc", default = "l2fc.csv", help = "Path to file containing l2fc values.")
parser$add_argument("--growth_pattern_col", default = "growth_pattern",
                    help = "Column containing growth pattern annotations.")
parser$add_argument("--l2fc_col", default = "l2fc", help = "Column containing log2 fold change values to be corrected.")
parser$add_argument("--sig_cols", default = "cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")

args = parser$parse_args()

# Read in files and set up any args
l2fc = data.table::fread(args$l2fc)
growth_pattern_col = args$growth_pattern_col
l2fc_col = args$l2fc_col
sig_cols = unlist(strsplit(args$sig_cols, ","))

# Detect if bio_rep column exists - MAYBE make this into a parameter?
if ("bio_rep" %in% colnames(l2fc)) {
  bio_rep_id_cols = c(sig_cols, "bio_rep")
} else {
  bio_rep_id_cols = sig_cols
  message("bio_rep column not detected. Assuming that there are NO biological replicates.")
}

# Drop cell set from grouping cols if is there
if ("cell_set" %in% bio_rep_id_cols) {
  message("Detecting cell_set column in bio_rep_id_cols. Dropping this column from the list.")
  bio_rep_id_cols = bio_rep_id_cols[bio_rep_id_cols != "cell_set"]
}

# Throw an error if a cell line does not have a growth annotation
if (any(unique(l2fc[[growth_pattern_col]]) %in% c(NA, "", " ", "NA"))) {
  stop("Detecting empty or missing values in the growth pattern column - ", growth_pattern_col)
}

# Correct l2fcs by regressing out cell line growth patterns
corrected_l2fc = l2fc |>
  dplyr::mutate(negcon_log2_norm_n = log2(control_median_normalized_n)) |>
  dplyr::group_split(dplyr::across(tidyselect::all_of(bio_rep_id_cols))) |>
  lapply(apply_bias_correction,
         raw_l2fc_col = l2fc_col,
         growth_pattern_col = growth_pattern_col,
         cell_set_col = "cell_set") |>
  dplyr::bind_rows() |>
  dplyr::ungroup()

# Check the column has been generated
if ("l2fc_uncorrected" %in% colnames(corrected_l2fc)) {
  if (!NA %in% unique(corrected_l2fc$l2fc_uncorrected)) {
    message("Success: l2fc_uncorrected column was generated.")
  } else {
    message("Warning: l2fc_uncorrected column was generated but with NA values.")
  }
} else {
  stop("Error: l2fc_uncorrected was not generated.")
}

# Write out
message("Writing corrected l2fc value to ", args$l2fc)
write.csv(corrected_l2fc, args$l2fc, row.names = FALSE, quote = FALSE)

# Check file creation
check_file_exists(args$l2fc)