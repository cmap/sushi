library(argparse)
library(magrittr)
library(tidyverse)

library(rhdf5)

source("./src/kitchen_utensils.R")
source("./src/compute_l2fc.R")
source("./src/synergy.R")


# TO DO: Added arguements ----

# Restructure L2FCs with single agent L2FCs ----
restructured_l2fc = restructure_l2fc(cps_l2fc = l2fc,
                                     cell_line_cols = cell_line_cols,
                                     sig_cols = sig_cols,
                                     singles_type = singles_type,
                                     combos_type = combos_type,
                                     l2fc_col = l2fc_col)

# Calculate synergies for combination treatments ----
trt_synergy = calculate_synergy(restructured_l2fc = restructured_l2fc,
                                l2fc_root = l2fc_col,
                                viability_cap = cap_for_viability)

# Calculate L2FC for DMSO ----
dmso_normalized_counts = normalized_counts %>%
  dplyr::filter(pert_type == "ctl_vehicle", pool_id != "CTLBC", pert_plate == "PCPS023")

dmso_dups = dmso_normalized_counts %>% dplyr::mutate(pert_type = "trt_cp") %>%
  dplyr::bind_rows(dmso_normalized_counts) %>%
  dplyr::select(-bio_rep) # will need to drop bio_rep column so that those are not collapsed

dmso_l2fc = compute_l2fc(normalized_counts = dmso_dups,
                         control_type = control_type,
                         sig_cols = c(id_cols, sig_cols),
                         ctrl_cols= c("pert_plate"),
                         count_col_name = 'logMFI_norm',
                         count_threshold = 0,
                         cell_line_cols= cell_line_cols)

# Loop to sample DMSO L2FC and write to an h5 ----
# Create an hf file
rhdf5::h5createFile("~/test_mock_synergy.h5")

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
  subset_combo = sample_dmso_l2fc(dmso_l2fc = subset,
                                  l2fc_col = "l2fc",
                                  grouping_cols = c("group_name"),
                                  n_samples = n_samples)
  
  # calculate synergy
  subset_synergy = calculate_synergy(restructured_l2fc = subset_combo,
                                     l2fc_root = "l2fc")
  
  # write to hdf5
  rhdf5::h5write(obj = subset_synergy$synergy, 
                 file = "~/Documents/cps_synergy_testing/test_mock_synergy.h5", 
                 name = i)
}

rhdf5::h5ls("~/test_mock_synergy.h5")
rhdf5::h5closeAll()
# by plate is resulting large dfs in memory

# Get pvalues ----
# Create group name and filter for only cls that were sampled
trt_synergy_scores_id = trt_synergy %>%
  tidyr::unite(all_of(grouping_cols), col = 'group_name', sep = "___", remove = FALSE, na.rm = FALSE)

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

# End ----

