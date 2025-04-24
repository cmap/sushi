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
parser$add_argument("--sig_cols", default = "cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "Column names that describe a unique profile")
parser$add_argument("--l2fc_col", default = "median_l2fc", help = "Name of column containing l2fc values")
parser$add_argument("--viab_cap", default = 1, help = "Maximum value for the viability")
parser$add_argument("--n_samples", default = 10000, help = "Number of samples used to create the null distribution")
parser$add_argument("-o", "--out", default = "", help = "Output path, defaults to working directory")

args = parser$parse_args()
#

combos_type = "trt_combo"
control_type = "ctl_vehicle"
dose_cols = c("pert1_dose", "pert2_dose")
singles_type = "trt_cp"

# Read in files and set up parameters ----
cps_l2fc = data.table::fread(args$l2fc, header = TRUE, sep = ",")
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
sig_cols = unlist(strsplit(args$sig_cols, ","))
#

# Create dmso l2fc ----
dmso_dups= dmso_entries %>% dplyr::mutate(pert_type = "trt_cp") %>%
  dplyr::bind_rows(dmso_entries) %>%
  dplyr::select(-bio_rep) # will need to drop bio_rep column so that those are not collapsed

dmso_l2fc = compute_l2fc(normalized_counts = dmso_dups,
                         control_type = control_type,
                         sig_cols = c(id_cols, sig_cols),
                         ctrl_cols= c("pert_plate"),
                         count_col_name = 'logMFI_norm',
                         count_threshold = 0,
                         cell_line_cols= cell_line_cols)
#

# Loop to sample and write to an hdf5 ----
# Create an hdf file
rhdf5::h5createFile("")

# Create unique group names for each plate cell line
# this is used for hdf5 hierarch
group_named_dmso_l2fc = dmso_l2fc %>% 
  tidyr::unite(all_of(grouping_cols), col = 'group_name', sep = "___", remove = TRUE, na.rm = FALSE)

# Loop through each unique name - large time constraint
unique_group_names = unique(group_named_dmso_l2fc$group_name)
for (i in unique_group_names) {
  print(i)
  
  subset = group_named_dmso_l2fc %>% dplyr::filter(group_name == i)
  
  # sample to n_samples and create pert1, pert2, combo
  resampled_l2fc = median_sample(x = subset$l2fc, n_samples = n_samples, size = 3, replace = FALSE, seed = 100)
  mock_values = data.table(pert1_l2fc = resampled_l2fc)
  mock_values[, pert2_l2fc := sample(pert1_l2fc, size = n_samples, replace = TRUE)]
  mock_values[, combo_l2fc := sample(pert1_l2fc, size = n_samples, replace = TRUE)]
  
  # calculate synergy
  subset_synergy = calculate_synergy(restructured_l2fc = mock_values, l2fc_root = "l2fc")
  
  
  # write to hdf5
  rhdf5::h5write(obj = subset_synergy$synergy, 
                 file = "", 
                 name = i)
}

rhdf5::h5ls("")
rhdf5::h5closeAll()
# by plate is resulting large dfs in memory

# Restructure L2FCs with single agent L2FCs ----
restructured_l2fc = restructure_l2fc(cps_l2fc = cps_l2fc,
                                     sig_cols = sig_cols,
                                     cell_line_cols = cell_line_cols,
                                     singles_type = singles_type,
                                     combos_type = combos_type,
                                     l2fc_col = l2fc_col,
                                     names_prefix = c("pert1", "pert2", "combo"),
                                     names_sep = "_")
#

# Get synergy scores ----
# Adds capped viability for all three l2fc columns
# Calculates hsa, bliss, and synergy
trt_synergy_scores = calculate_synergy(restructured_l2fc = restructured_l2fc,
                                       l2fc_root = l2fc_col,
                                       viability_cap = cap_for_viability)
#

# Pull out pvalues ----
# Some celll lines may not be present. hdf5 creation was terminated early
# Open hdf5 connection
test_mock_synergy = rhdf5::H5Fopen("")
temp = rhdf5::h5ls("")

# Create group name and filter for only cls that were sampled
trt_synergy_scores_id = trt_synergy_scores %>%
  tidyr::unite(all_of(grouping_cols), col = 'group_name', sep = "___", remove = FALSE, na.rm = FALSE) %>%
  dplyr::filter(group_name %in% temp$name)

# Is there anyway to do this in dplyr instead of as an apply?
# This doesnt take that long
test = apply(trt_synergy_scores_id[,c('group_name','synergy')], 1, 
             function(y) smooth_pvalue(get_cdf_value(y['group_name'] ,y['synergy'], h5_file = test_mock_synergy),
                                       n_samples = n_samples))
rhdf5::h5closeAll()

trt_synergy_scores_id %<>% 
  dplyr::mutate(cdf_val = as.numeric(test)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(p_val.emp = 2 * min(1 - cdf_val, cdf_val)) %>% 
  dplyr::ungroup()

# Write out file ----
#collapsed_l2fc_outpath = paste(args$out, 'collapsed_l2fc.csv', sep='/')
#print(paste0('Writing out collapsed l2fc file to ', collapsed_l2fc_outpath))
#write.csv(x= collapsed_l2fc, file= collapsed_l2fc_outpath, row.names= FALSE, quote= FALSE)
#

# End ----