library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
library(rhdf5)

source("utils/kitchen_utensils.R")
source("compute_l2fc/compute_l2fc_functions.R")
source("synergy/synergy_functions.R")

# Shell script argument parser ----
parser = ArgumentParser()

parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE, help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose", help = "Print little output")
parser$add_argument("--wkdir", default = getwd(), help = "Working directory")
parser$add_argument("--normalized_counts", default = "normalized_counts.csv",
                    help = "Path to file containing normalized counts")
parser$add_argument("--l2fc", default = "l2fc.csv", help = "Path to file containing l2fc values")
# Used to create mock DMSO l2fcs and to restructure l2fcs
parser$add_argument("--cell_line_cols", default = "pool_id,depmap_id,lua",
                    help = "Column names that describe a cell line")
# Used to create mock DMSO l2fcs
parser$add_argument("--ctrl_cols", default = "cell_set,day",
                    help = "Column names that describe the control conditions")
# Used to create mock DMSO l2fcs and to restructure l2fcs
parser$add_argument("--sig_cols", default = "cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "Column names that describe a unique profile")
# Used to restructure l2fcs
parser$add_argument("--combination_col", default = "is_combination", help = "Column name used to identify combinations")
# Used to create mock DMSO l2fcs
parser$add_argument("--count_col_name", default = "log2_normalized_n",
                    help = "Column name containing counts to be used to calculate l2fcs")
# Used to create mock DMSO l2fcs
parser$add_argument("--count_threshold", default = 40, help = "Low counts threshold")
# Used to restructure l2fcs and to calculate synergy
parser$add_argument("--l2fc_col", default = "median_l2fc", help = "Column name containing l2fc values")
# Used to create mock DMSO l2fcs
parser$add_argument("--n_samples", default = 10000, help = "Size of the resampling")
# Used to create mock DMSO l2fcs
parser$add_argument("--negcon_type", default = "ctl_vehicle",
                    help = "String in pert_type that identifies the negative control")
# Used to calculate synergy
parser$add_argument("--viab_cap", default = 1, help = "Maximum value for the viability")
parser$add_argument("-o", "--out", default = "", help = "Output path, defaults to working directory")

args = parser$parse_args()

# Read in files and set up parameters ----
normalized_counts = data.table::fread(args$normalized_counts, sep = ",", header = TRUE)
cps_l2fc = data.table::fread(args$l2fc, sep = ",", header = TRUE)

args$cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
args$ctrl_cols = unlist(strsplit(args$ctrl_cols, ","))
args$sig_cols = unlist(strsplit(args$sig_cols, ","))
args$count_threshold = as.numeric(args$count_threshold)
args$n_samples = as.numeric(args$n_samples)
args$viab_cap = as.numeric(args$viab_cap)

# Create DMSO l2fc ----
dmso_norm_ctrl = normalized_counts[pert_type == args$negcon_type, ]
dmso_norm_trt = dmso_norm_ctrl
dmso_norm_trt$pert_type = "trt_cp"
dmso_norm_input = data.table::rbindlist(list(dmso_norm_trt, dmso_norm_ctrl))

dmso_l2fc = compute_l2fc(normalized_counts = dmso_norm_input,
                         cell_line_cols = args$cell_line_cols,
                         control_type = args$negcon_type,
                         count_col_name = args$count_col_name,
                         count_threshold = args$count_threshold,
                         ctrl_cols = args$ctrl_cols,
                         sig_cols = args$sig_cols)

# Loop to resample and write to an h5 ----
# Create an hdf file
rhdf5::h5createFile(file.path(args$out, "dmso_synergy.h5"))

# Create unique group names for each plate + cell line for hdf5 hierarchy
dmso_l2fc = data.table::as.data.table(dmso_l2fc)
group_name_cols = unique(c(args$cell_line_cols, args$ctrl_cols)) # also used later for mapping
dmso_l2fc[, group_name := do.call(paste, c(.SD, sep = "___")), .SDcols = group_name_cols]

# Loop through each unique name
unique_group_names = unique(dmso_l2fc$group_name)
for (i in unique_group_names) {
  subset = dmso_l2fc[group_name == i, ]

  # Resample to n_samples and create pert1, pert2, combo
  resampled_l2fc = median_resample(x = subset$l2fc, n_samples = args$n_samples, size = 3, replace = FALSE, seed = 2)
  mock_values = data.table(pert1_l2fc = resampled_l2fc)
  mock_values[, pert2_l2fc := sample(pert1_l2fc, size = args$n_samples, replace = TRUE)]
  mock_values[, combo_l2fc := sample(pert1_l2fc, size = args$n_samples, replace = TRUE)]

  # Calculate synergy - function occurs in place
  calculate_synergy(restructured_l2fc = mock_values,
                    l2fc_cols = colnames(mock_values))
  # Write to hdf5
  rhdf5::h5write(obj = mock_values$synergy,
                 file = file.path(args$out, "dmso_synergy.h5"),
                 name = i)
}

print("Populated dmso_synergy.h5.")
rhdf5::h5closeAll()

# Restructure L2FCs with single agent L2FCs and calculate synergy ----
# Create vector of new column names
names_prefix = c("pert", "pert2", "combo")
new_names = paste(names_prefix, args$l2fc_col, sep = "_")

# Pull out columns that describe each perturbation
ignore_cols = c("pert_type", "pert_plate") # common columns with "pert" to be ignored
pert1_cols = args$sig_cols[grepl(paste0(names_prefix[1], "_"), args$sig_cols) & !args$sig_cols %in% ignore_cols]
pert2_cols = args$sig_cols[grepl(paste0(names_prefix[2], "_"), args$sig_cols) & !args$sig_cols %in% ignore_cols]

# Pull columns to join on - these are sig_cols that don't describe the perturbations
filt_sig_cols = args$sig_cols[!args$sig_cols %in% c(pert1_cols, pert2_cols)]

trt_synergy = restructure_l2fc(cps_l2fc = cps_l2fc,
                               join_cols = unique(c(args$cell_line_cols, filt_sig_cols)),
                               pert_cols_list = list(pert1_cols, pert2_cols),
                               l2fc_col = args$l2fc_col,
                               new_l2fc_col_names = new_names)

# Add synergy scores to data.table in-place
calculate_synergy(restructured_l2fc = trt_synergy,
                  l2fc_cols = new_names,
                  viab_cap = as.numeric(args$viab_cap))
print("Calculated synergies.")

# Pull out pvalues ----
# Open hdf5 connection
dmso_synergy = rhdf5::H5Fopen(file.path(args$out, "dmso_synergy.h5"))
existing_groups = rhdf5::h5ls(file.path(args$out, "dmso_synergy.h5"))

# Create group name and filter for only cls that were sampled
trt_synergy[, group_name := do.call(paste, c(.SD, sep = "___")), .SDcols = group_name_cols]
# Filter out names not in hf file
trt_synergy = trt_synergy[group_name %in% existing_groups$name, ]

# Mapply vectorize
trt_synergy[, p_val_emp := mapply(get_pvalue, group_name, synergy,
                                  MoreArgs = list(h5_file = dmso_synergy, n_samples = args$n_samples))]
rhdf5::h5closeAll()
print("Added empirical pvalues to synergies.")

# Write out file ----
outpath = file.path(args$out, "synergy_scores.csv")
print(paste0("Writing out synergy file to ", outpath))
write.csv(x = trt_synergy, file = outpath, row.names = FALSE)

# End ----