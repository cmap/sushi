options(cli.unicode = FALSE)
library(tidyverse)
library(data.table)
library(magrittr)
library(data.table)
library (PRROC)

# BY ID_COLS (PCR_PLATE, PCR_WELL) ----------

# Compute skew
compute_skew <- function(df, group_cols, metric) {
  # Compute skew based on AUC
  df %>%
    filter(expected_read) %>%
    group_by(across(all_of(group_cols))) %>%
    arrange(desc(.data[[metric]])) %>%  # Sort by metric in descending order
    mutate(
      rank_fraction = row_number() / n(),  # X-axis: rank fraction
      cumulative_fraction_reads = cumsum(.data[[metric]]) / sum(.data[[metric]], na.rm = TRUE)  # Y-axis: cumulative fraction
    ) %>%
    summarise(
      skew = if (n() > 1) {
        # Compute AUC using trapezoidal rule
        sum((rank_fraction[-1] - rank_fraction[-n()]) *
              (cumulative_fraction_reads[-1] + cumulative_fraction_reads[-n()]) / 2, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

compute_expected_lines <- function(cell_set_meta) {
  # Get number of expected cell lines for each cell_set
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

compute_read_stats <- function(annotated_counts, cell_set_meta, group_cols = c("pcr_plate","pcr_well"),
                               metric = "n") {
  # Compute expected lines from cell_set_meta
  expected_lines <- compute_expected_lines(cell_set_meta)
  # Convert to a dataframe
  expected_lines <- as.data.table(expected_lines)

  result <- annotated_counts %>%
    left_join(expected_lines, by = "cell_set") %>% # Add n_expected_lines from lookup
    group_by(across(all_of(group_cols))) %>%
    summarise(
      # Total reads
      n_total_reads = sum(.data[[metric]], na.rm = TRUE),
      # Reads mapping to cell lines
      n_cl_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      # Fraction of reads mapping to cell lines
      fraction_expected_reads = n_cl_reads / n_total_reads,
      # Reads mapping to control barcodes
      n_cb_reads = sum(!is.na(cb_name) & cb_name != "", na.rm = TRUE),
      # Spearman correlation of control barcodes
      cb_spearman = if (sum(!is.na(cb_name) & !is.na(cb_log2_dose)) > 1) cor(
        cb_log2_dose[!is.na(cb_name) & !is.na(cb_log2_dose)],
        .data[[metric]][!is.na(cb_name) & !is.na(cb_log2_dose)],
        method = "spearman",
        use = "complete.obs"
      ) else NA_real_,
      # Number of cell lines with coverage below 50 and 90 reads
      n_cl_below_50 = sum(.data[[metric]] < 50, na.rm = TRUE),
      n_cl_below_95 = sum(.data[[metric]] < 90, na.rm = TRUE),
      # Number of cell lines with coverage above 40 reads
      n_lines_recovered = sum(.data[[metric]] >= 40 & (is.na(cb_name) | cb_name == ""), na.rm = TRUE),
      # Number of expected lines based on metadata
      n_expected_lines = max(n_expected_lines, na.rm = TRUE), # Bring forward from join
      # Fraction of cell lines with coverage above 40 reads
      fraction_cl_recovered = n_lines_recovered / max(n_expected_lines, na.rm = TRUE),
      .groups = "drop"
    )
}

# CELL LINE BY PLATE (PCR_PLATE, CELL LINE) ----------

# Compute error rate
compute_error_rate <- function(df, metric = 'log2_normalized_n', group_cols = c("depmap_id", "pcr_plate"),
                               negcon = "ctl_vehicle", poscon = "trt_poscon") {
  df %>%
    filter(
      pert_type %in% c(negcon, poscon),
      is.finite(.data[[metric]]),
      !is.na(pool_id)
    ) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      error_rate = {
        roc_data <- PRROC::roc.curve(
          scores.class0 = .data[[metric]],
          weights.class0 = pert_type == negcon,
          curve = TRUE
        )
        min(roc_data$curve[, 1] + 1 - roc_data$curve[, 2]) / 2
      },
      .groups = "drop"
    )
}


# Compute vehicle and positive control median and MAD for both raw and normalized data
compute_ctl_medians_and_mad <- function(df, group_cols = c("depmap_id", "pcr_plate"),
                                    negcon = "ctl_vehicle", poscon = "trt_poscon") {
  # Group and compute medians/MADs
  df %>%
    filter(pert_type %in% c(negcon, poscon)) %>%
    group_by(across(all_of(c(group_cols, "pert_type")))) %>%
    summarise(
      median_normalized = median(log2_normalized_n, na.rm = TRUE),
      mad_normalized = mad(log2_normalized_n, na.rm = TRUE),
      median_raw = median(log2_n, na.rm = TRUE),
      mad_raw = mad(log2_n, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = pert_type,
      values_from = c(median_normalized, mad_normalized, median_raw, mad_raw),
      names_sep = "_"
    )
}

# Compute log fold change
compute_control_lfc <- function(df, negcon = "ctl_vehicle", poscon = "trt_poscon") {
  df %>%
    mutate(
      lfc_normalized = .data[[paste0("median_normalized_", poscon)]] -
                       .data[[paste0("median_normalized_", negcon)]],
      lfc_raw = .data[[paste0("median_raw_", poscon)]] -
                .data[[paste0("median_raw_", negcon)]]
    )
}

compute_false_sensitivity <- function(df, negcon_type = "ctl_vehicle") {
  # Define column to use
  col <- paste0("mad_normalized_", negcon_type)
  print(paste0("Computing false sensitivity probability using column: ", col))

  # Group by the specified columns and calculate the probability
  result <- df %>%
    mutate(
      false_sensitivity_probability = pnorm(
        -2,
        sd = .data[[col]] * sqrt((1/32 + 1/3) * pi/2)
      )
    )

  return(result)
}

compute_cl_fractions <- function(df, metric = "n", grouping_cols = c("pcr_plate", "depmap_id")) {
  df %>%
    group_by(across(all_of(grouping_cols))) %>%
    mutate(
      total_reads = sum(.data[[metric]], na.rm = TRUE),  # Total reads per group
      fraction_of_reads = .data[[metric]] / total_reads  # Fraction of reads for each entry
    ) %>%
    ungroup() %>%
    select(all_of(grouping_cols), !!metric)  # Retain relevant columns
}