options(cli.unicode = FALSE)
library(tidyverse)
library(data.table)
library(magrittr)
library(data.table)
library (PRROC)

# BY ID_COLS (PCR_PLATE, PCR_WELL) ----------

# Compute skew
compute_skew <- function(df, group_cols = c("pcr_plate","pcr_well"), metric = "n") {
  df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    arrange(desc(.data[[metric]])) %>%  # Sort by metric in descending order
    mutate(
      rank_fraction = row_number() / n(),  # Calculate rank fraction of each cell line
      cumulative_fraction_reads = cumsum(.data[[metric]]) / sum(.data[[metric]], na.rm = TRUE) # Calculate cumulative fraction of reads
    ) %>%
    summarise(
      skew = if (n() > 1) {
        # Compute AUC
        sum((rank_fraction[-1] - rank_fraction[-n()]) *
              (cumulative_fraction_reads[-1] + cumulative_fraction_reads[-n()]) / 2, na.rm = TRUE)
      }
    ) %>%
    dplyr::ungroup()
}

compute_expected_lines <- function(cell_set_meta) {
  # Get number of expected cell lines for each cell_set
  cell_set_meta %>%
    dplyr::group_by(depmap_id) %>%
    dplyr::filter(n() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cell_set) %>%
    dplyr::summarise(
      n_expected_lines = dplyr::n_distinct(depmap_id), # Count unique depmap_id for each cell_set
    ) %>%
    dplyr::ungroup()
}

compute_read_stats <- function(annotated_counts, cell_set_meta, group_cols = c("pcr_plate","pcr_well"),
                               metric = "n") {
  # Compute expected lines from cell_set_meta
  expected_lines <- compute_expected_lines(cell_set_meta)

  result <- annotated_counts %>%
    dplyr::left_join(expected_lines, by = "cell_set") %>% # Add n_expected_lines from lookup
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise(
      # Total reads
      n_total_reads = dplyr::sum(.data[[metric]], na.rm = TRUE),
      # Reads mapping to cell lines
      n_cl_reads = dplyr::sum(.data[[metric]][expected_read], na.rm = TRUE),
      # Fraction of reads mapping to cell lines
      fraction_expected_reads = n_cl_reads / n_total_reads,
      # Reads mapping to control barcodes
      n_cb_reads = dplyr::sum(!is.na(cb_name) & cb_name != "", na.rm = TRUE),
      # Spearman correlation of control barcodes
      cb_spearman = if (dplyr::sum(!is.na(cb_name) & !is.finite(cb_log2_dose)) > 1) cor(
        cb_log2_dose[!is.na(cb_name) & !is.finite(cb_log2_dose)],
        .data[[metric]][!is.na(cb_name) & !is.finite(cb_log2_dose)],
        method = "spearman",
        use = "complete.obs"
      ),
      # Number of cell lines with coverage below 50 and 90 reads
      n_cl_below_50 = dplyr::sum(.data[[metric]] < 50, na.rm = TRUE),
      n_cl_below_95 = dplyr::sum(.data[[metric]] < 90, na.rm = TRUE),
      # Number of cell lines with coverage above 40 reads
      n_lines_recovered = dplyr::sum(.data[[metric]] >= 40 & (is.na(cb_name) | cb_name == ""), na.rm = TRUE),
      # Number of expected lines based on metadata
      n_expected_lines = max(n_expected_lines, na.rm = TRUE), # Bring forward from join
      # Fraction of cell lines with coverage above 40 reads
      fraction_cl_recovered = n_lines_recovered / max(n_expected_lines, na.rm = TRUE),
    ) %>%
    dplyr::ungroup()
}

# CELL LINE BY PLATE (PCR_PLATE, CELL LINE) ----------

# Compute error rate
compute_error_rate <- function(df, metric = 'log2_normalized_n', group_cols = c("depmap_id", "pcr_plate"),
                               negcon = "ctl_vehicle", poscon = "trt_poscon") {
  paste0("Computing error rate using ", negcon, " and ", poscon, ".....")
  paste0("Grouping by ", paste0(group_cols, collapse = ","), ".....")
  df %>%
    dplyr::filter(
      pert_type %in% c(negcon, poscon),
      is.finite(.data[[metric]]),
      !is.na(pool_id)
    ) %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise(
      error_rate = {
        roc_data <- PRROC::roc.curve(
          scores.class0 = .data[[metric]],
          weights.class0 = pert_type == negcon,
          curve = TRUE
        )
        min(roc_data$curve[, 1] + 1 - roc_data$curve[, 2]) / 2
      },
    ) %>%
    dplyr::ungroup()
}


# Compute vehicle and positive control median and MAD for both raw and normalized data
compute_ctl_medians_and_mad <- function(df, group_cols = c("depmap_id", "pcr_plate"),
                                    negcon = "ctl_vehicle", poscon = "trt_poscon") {
  paste0("Adding control median and MAD values for ", negcon, " and ", poscon, ".....")
  paste0("Computing falses sensitivity probability for ", negcon, ".....")
  # Group and compute medians/MADs
  df %>%
    dplyr::filter(pert_type %in% c(negcon, poscon)) %>%
    dplyr::group_by(across(all_of(c(group_cols, "pert_type")))) %>%
    dplyr::summarise(
      median_normalized = median(log2_normalized_n, na.rm = TRUE),
      n_replicates = n(),
      mad_normalized = mad(log2_normalized_n, na.rm = TRUE),
      median_raw = median(log2_n, na.rm = TRUE),
      mad_raw = mad(log2_n, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    pivot_wider(
      names_from = pert_type,
      values_from = c(median_normalized, mad_normalized, median_raw, mad_raw, n_replicates),
      names_sep = "_"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      false_sensitivity_probability_50 = pnorm(
        -1,
        sd = .data[[paste0('mad_normalized_',negcon)]] * sqrt((1/.data[[paste0('n_replicates_',negcon)]] * pi/2) + 1)
      ),
      false_sensitivity_probability_25 = pnorm(
        -2,
        sd = .data[[paste0('mad_normalized_',negcon)]] * sqrt((1/.data[[paste0('n_replicates_',negcon)]] * pi/2) + 1)
      )
    )
}

# Compute log fold change
compute_control_lfc <- function(df, negcon = "ctl_vehicle", poscon = "trt_poscon") {
  paste0("Computing log fold change for ", negcon, " and ", poscon, ".....")
  df %>%
    dplyr::mutate(
      lfc_normalized = .data[[paste0("median_normalized_", poscon)]] -
                       .data[[paste0("median_normalized_", negcon)]],
      lfc_raw = .data[[paste0("median_raw_", poscon)]] -
                .data[[paste0("median_raw_", negcon)]]
    )
}

compute_cl_fractions <- function(df, metric = "n", grouping_cols = c("pcr_plate", "depmap_id")) {
  paste0("Computing cell line fractions for ", metric, ".....")
  df %>%
    dplyr::group_by(across(all_of(grouping_cols))) %>%
    dplyr::mutate(
      total_reads = sum(.data[[metric]], na.rm = TRUE),  # Total reads per group
      fraction_of_reads = .data[[metric]] / total_reads  # Fraction of reads for each entry
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(grouping_cols), !!metric)  # Retain relevant columns
}