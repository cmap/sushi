library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
library(rhdf5)

source("./src/kitchen_utensils.R")
source("./src/compute_l2fc.R")
source("./src/synergy.R")

# Shell script argument parser ----
parser = ArgumentParser()

parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE, help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose", help = "Print little output")
parser$add_argument("--wkdir", default = getwd(), help = "Working directory")
parser$add_argument("--l2fc", default = "l2fc.csv", help = "Path to file containing l2fc values")
parser$add_argument("--cell_line_cols", default = "pool_id,depmap_id,lua",
                    help = "Column names that describe a cell line")
parser$add_argument("--ctrl_cols", default="cell_set,day", 
                    help = "columns used to collapse controls to generate l2fc")
parser$add_argument("--control_type", default="ctl_vehicle", help="pert_type to use as control")
parser$add_argument("--count_col_name", default="log2_normalized_n", 
                    help = "column containing counts with which to calculate l2fc")
parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
parser$add_argument("--sig_cols", default = "cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "Column names that describe a unique profile")
parser$add_argument("--l2fc_col", default = "median_l2fc", help = "Name of column containing l2fc values")
parser$add_argument("--viab_cap", default = 1, help = "Maximum value for the viability")
parser$add_argument("--n_samples", default = 10000, help = "Number of samples used to create the null distribution")
parser$add_argument("-o", "--out", default = "", help = "Output path, defaults to working directory")

args = parser$parse_args()
#

# Read in files and set up parameters ----
normalized_counts = data.table::fread(args$normalized_counts, sep = ",", header = TRUE)
cps_l2fc = data.table::fread(args$l2fc, sep = ",", header = TRUE)

args$cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
args$count_threshold = as.numeric(args$count_threshold)
args$ctrl_cols = unlist(strsplit(args$ctrl_cols, ","))
args$n_samples = as.numeric(args$n_samples)
args$sig_cols = unlist(strsplit(args$sig_cols, ","))
args$viab_cap = as.numeric(args$viab_cap)
#

# Create DMSO l2fc ----
dmso_norm_ctrl = normalized_counts[pert_type == args$control_type, ]
dmso_norm_trt = dmso_norm_ctrl
dmso_norm_trt$pert_type = "trt_cp"
dmso_norm_input = data.table::rbindlist(list(dmso_norm_trt, dmso_norm_ctrl))

dmso_l2fc = compute_l2fc(normalized_counts = dmso_norm_input,
                         cell_line_cols = args$cell_line_cols,
                         control_type = args$control_type,
                         count_col_name = args$count_col_name,
                         count_threshold = args$count_threshold,
                         ctrl_cols = args$ctrl_cols,
                         sig_cols = args$sig_cols)
#

# Loop to resample and write to an hdf5 ----
# Create an hdf file
rhdf5::h5createFile(file.path(args$out, "mock_dmso_synergy.h5"))

# Create unique group names for each plate cell line
# this is used for hdf5 hierarchy
dmso_l2fc = data.table::as.data.table(dmso_l2fc)
dmso_l2fc[, group_name := do.call(paste, c(.SD, sep = "___")), .SDcols = c(args$cell_line_cols, args$ctrl_cols)]

# Loop through each unique name - large time constraint
unique_group_names = unique(dmso_l2fc$group_name)
for (i in unique_group_names) {
  print(i)
  subset = dmso_l2fc[group_name == i, ]

  # Resample to n_samples and create pert1, pert2, combo
  resampled_l2fc = median_sample(x = subset$l2fc, n_samples = args$n_samples, size = 4,
                                 replace = FALSE, seed = 100)
  mock_values = data.table(pert1_l2fc = resampled_l2fc)
  mock_values[, pert2_l2fc := sample(pert1_l2fc, size = args$n_samples, replace = TRUE)]
  mock_values[, combo_l2fc := sample(pert1_l2fc, size = args$n_samples, replace = TRUE)]

  # Calculate synergy
  subset_synergy = calculate_synergy(restructured_l2fc = mock_values,
                                     l2fc_cols = c("pert1_l2fc", "pert2_l2fc", "combo_l2fc"))
  # Write to hdf5
  rhdf5::h5write(obj = subset_synergy$synergy,
                 file = file.path(args$out, "mock_dmso_synergy.h5"),
                 name = i)
}

rhdf5::h5closeAll()
# by plate is resulting large dfs in memory

# Restructure L2FCs with single agent L2FCs and calculate synergy ----
# Create vector of new column names
names_prefix = c("pert", "pert2", "combo")
new_names = paste(names_prefix, args$l2fc_col, sep = "_")

# Pull out columns that describe each pert
ignore_cols = c("pert_type", "pert_plate") # common pert columns that could be extracted
pert1_cols = args$sig_cols[grepl(paste0(names_prefix[1], "_"), args$sig_cols) & !args$sig_cols %in% ignore_cols]
pert2_cols = args$sig_cols[grepl(paste0(names_prefix[2], "_"), args$sig_cols) & !args$sig_cols %in% ignore_cols]

# Pull columns to join on
# These are cell line columns and additional columns from sig_cols that don't describe the pertubation
filt_sig_cols = args$sig_cols[!args$sig_cols %in% c(pert1_cols, pert2_cols, "pert_type")]

restructured_l2fc = restructure_l2fc(cps_l2fc = cps_l2fc,
                                     join_cols = unique(c(args$cell_line_cols, filt_sig_cols)),
                                     pert_cols_list = list(pert1_cols, pert2_cols),
                                     l2fc_col = args$l2fc_col,
                                     single_type = "trt_cp",
                                     combo_type = "trt_combo",
                                     new_col_names = new_names)

trt_synergy_scores = calculate_synergy(restructured_l2fc = restructured_l2fc,
                                       l2fc_cols = new_names,
                                       viab_cap = as.numeric(args$viab_cap))
#

# Pull out pvalues ----
# Open hdf5 connection
mock_dmso_synergy = rhdf5::H5Fopen(file.path(args$out, "mock_dmso_synergy.h5"))
existing_groups = rhdf5::h5ls(file.path(args$out, "mock_dmso_synergy.h5"))

# Create group name and filter for only cls that were sampled
trt_synergy_scores[, group_name := do.call(paste, c(.SD, sep = "___")),
                   .SDcols = c(args$cell_line_cols, args$ctrl_cols)]
# Filter out names not in hf file
trt_synergy_scores = trt_synergy_scores[group_name %in% existing_groups$name, ]

# Mapply vectorize
trt_synergy_scores[, p_val_emp := mapply(get_pvalue, group_name, synergy,
                                         MoreArgs = list(h5_file = mock_dmso_synergy,
                                                         n_samples = args$n_samples))]
rhdf5::h5closeAll()
#

# Write out file ----
outpath = file.path(args$out, "synergy_scores.csv")
print(paste0("Writing out synergy file to ", outpath))
write.csv(x = trt_synergy_scores, file = outpath, row.names = FALSE)
#

# End ----