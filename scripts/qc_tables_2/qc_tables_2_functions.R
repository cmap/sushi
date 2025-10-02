options(cli.unicode = FALSE)
library(tidyverse)
library(data.table)
library(magrittr)
library(PRROC)
library(dplyr)


# Pool treatment replicate concordance ----------

compute_pool_delta_df <- function(l2fc, threshold = 2.0) {
  # Define the grouping columns for the window function and the threshold
  grouping_cols <- c('day', 'pert_name', 'pert_dose', 'pert_type', 'pool_id', 'pert_plate', 'cell_set')
  DELTA_THRESHOLD <- threshold

  # NEW: Calculate the number of cell lines with low fold change for each pool.
  # This must be done before the main data is summarized.
  killed_df <- l2fc %>%
    mutate(fc = 2^l2fc) %>%
    group_by(pert_name, pert_dose) %>%
    summarise(
      # Count distinct cell lines where fc < 0.25
      n_lines_killed = n_distinct(depmap_id[fc < 0.25 & !is.na(fc)]),
      .groups = 'drop'
    )

  # Calculate the median l2fc and fc for each biological replicate pool
  pool_median_df <- l2fc %>%
    # Create the 'fc' column (2^l2fc)
    mutate(fc = 2^l2fc) %>%
    # Group by all replicate-specific columns
    group_by(across(all_of(grouping_cols)), bio_rep) %>%
    # Calculate the median for each pool
    summarise(
      pool_l2fc = median(l2fc, na.rm = TRUE),
      pool_fc = median(fc, na.rm = TRUE),
      .groups = 'drop' # Ungroup after summarising
    )

  # Calculate group-wide medians and deltas from those medians
  delta_df <- pool_median_df %>%
    # Group by the treatment-level columns (excluding bio_rep)
    group_by(across(all_of(grouping_cols))) %>%
    # The mutate() on a grouped data frame acts like a window function
    mutate(
      # Calculate and store the median across biological replicates for each group
      group_median_l2fc = median(pool_l2fc, na.rm = TRUE),
      group_median_fc = median(pool_fc, na.rm = TRUE)
    ) %>%
    ungroup() %>% # Ungroup to perform calculations on the whole data frame again
    mutate(
      # Use the new median column to calculate the delta
      delta_from_median_l2fc = pool_l2fc - group_median_l2fc,
      delta_from_median_fc = pool_fc - group_median_fc,
      # Flag rows where the delta exceeds the threshold
      is_outlier_l2fc = delta_from_median_l2fc > DELTA_THRESHOLD,
      is_outlier_fc = delta_from_median_fc > DELTA_THRESHOLD
    ) %>%
    # NEW: Join the fraction data back into the main results
    left_join(killed_df, by = c("pert_name", "pert_dose")) %>%
    # Filter out rows where pool_id is the string "NA"
    filter(pool_id != "NA") %>%
    select(
      pert_dose, pert_name, bio_rep, pert_plate, cell_set,
      delta_from_median_fc, delta_from_median_l2fc, is_outlier_fc, pool_id,
      group_median_fc, group_median_l2fc,
      n_lines_killed
    )

  return(delta_df)
}
