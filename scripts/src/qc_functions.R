library(tidyverse)
library(data.table)
library(magrittr)
library(data.table)

# BY ID_COLS (PCR_PLATE, PCR_WELL) ----------

# Helper functions ----------

# Compute skew
compute_skew <- function(df, group_cols, metric) {
  df %>%
    filter(expected_read) %>%
    group_by(across(all_of(group_cols))) %>%
    arrange(depmap_id) %>%
    mutate(
      rank_fraction = row_number() / n(),
      cumulative_fraction_reads = cumsum(.data[[metric]]) / sum(.data[[metric]], na.rm = TRUE)
    ) %>%
    summarise(
      auc = sum((rank_fraction[-1] - rank_fraction[-n()]) *
                  (cumulative_fraction_reads[-1] + cumulative_fraction_reads[-n()]) / 2, na.rm = TRUE),
      skew = 1 - auc,
      .groups = "drop"
    )
}

compute_expected_lines <- function(cell_set_meta) {
  cell_set_meta %>%
    group_by(depmap_id) %>%
    filter(n() == 1) %>%
    ungroup() %>%
    group_by(cell_set) %>%
    summarise(
      n_expected_lines = n_distinct(depmap_id), # Count unique depmap_id for each cell_set
      .groups = "drop"
    )
}

# Main function ----------
id_cols_table <- function(annotated_counts, cell_set_meta, group_cols, metric) {
  # Compute expected lines from cell_set_meta
  expected_lines <- compute_expected_lines(cell_set_meta)

  # Main annotated counts summarization
  result <- annotated_counts %>%
    left_join(expected_lines, by = "cell_set") %>% # Add n_expected_lines from lookup
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_total_reads = sum(.data[[metric]], na.rm = TRUE),
      n_cl_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      fraction_expected_reads = n_cl_reads / n_total_reads,
      n_cb_reads = sum(!is.na(cb_name) & cb_name != "", na.rm = TRUE),
      cb_spearman = if (sum(!is.na(cb_name) & !is.na(cb_log2_dose)) > 1) cor(
        cb_log2_dose[!is.na(cb_name) & !is.na(cb_log2_dose)],
        .data[[metric]][!is.na(cb_name) & !is.na(cb_log2_dose)],
        method = "spearman",
        use = "complete.obs"
      ) else NA_real_,
      n_cl_below_50 = sum(.data[[metric]] < 50, na.rm = TRUE),
      n_cl_below_95 = sum(.data[[metric]] < 90, na.rm = TRUE),
      n_lines_recovered = sum(.data[[metric]] >= 40 & (is.na(cb_name) | cb_name == ""), na.rm = TRUE),
      n_expected_lines = max(n_expected_lines, na.rm = TRUE), # Bring forward from join
      fraction_cl_recovered = n_lines_recovered / max(n_expected_lines, na.rm = TRUE),
      .groups = "drop"
    )

  # Compute skew
  skew <- compute_skew(annotated_counts, group_cols, metric)

  # Merge skew into result
  result <- result %>%
    left_join(skew, by = group_cols)

  return(result)
}


# Testing ----------

filtered_counts <- data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/filtered_counts.csv", header= TRUE, sep = ',')
sample_meta <- data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/sample_meta.csv", header= TRUE, sep = ',')
cell_set_meta <- data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/cell_set_and_pool_meta.csv", header= TRUE, sep = ',')
annotated_counts <- data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/annotated_counts.csv", header= TRUE, sep = ',')
normalized_counts <- data.table::fread("/Users/jdavis/Downloads/MTS_SEQ003/normalized_counts.csv", header= TRUE, sep = ',')

table_out <- id_cols_table(annotated_counts = annotated_counts, cell_set_meta = cell_set_meta, group_cols = c("pcr_plate", "pcr_well", "cell_set"), metric = "n")