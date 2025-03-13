options(cli.unicode = FALSE)
library(tidyverse)
library(data.table)
library(magrittr)
library (PRROC)

# BY ID_COLS (PCR_PLATE, PCR_WELL) ----------

#' Compute Skew
#'
#' This function computes the skew, which measures the cumulative fraction of barcode read counts taken up by each cell
#' line. It is computed as the auc of that CDF function and has a range of (0.5,1).
#' A lower skew indicates a more even distribution of reads across cell lines.
#'
#' @param df A data frame containing the data to compute skew, generally annotated_counts.
#' @param group_cols A character vector specifying the column names to group by (default: `c("pcr_plate", "pcr_well")`).
#' @param metric A string indicating the column name of the metric to use for calculations (default: `"n"`).
#'
#' @return A data frame with one row per group and a column `skew` containing the computed skew (auc).
#'
#' @import dplyr
compute_skew <- function(df, group_cols = c("pcr_plate","pcr_well"), metric = "n") {
  result <- df %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    arrange(desc(.data[[metric]])) %>%  # Sort by metric in descending order
    mutate(
      rank_fraction = row_number() / n(),  # Calculate rank fraction of each cell line
      cumulative_fraction_reads = cumsum(.data[[metric]]) / sum(.data[[metric]], na.rm = TRUE) # Calculate cumulative fraction of reads
    ) %>%
    summarise(
      skew = if (n() > 1) {
        # Compute auc
        sum((rank_fraction[-1] - rank_fraction[-n()]) *
              (cumulative_fraction_reads[-1] + cumulative_fraction_reads[-n()]) / 2, na.rm = TRUE)
      }
    ) %>%
    dplyr::ungroup()
  return(result)
}

#' Compute expected cell lines
#'
#' This function calculates the number of unique expected cell lines for each cell_set.
#'
#' @param cell_set_meta A data frame containing metadata for cell sets.
#' @param cell_line_cols A character vector specifying the column names that define a unique cell line.
#'
#' @return A data frame with one row per cell_set and a column `n_expected_lines` indicating the number of unique expected cell lines in each set.
#'
#' @import dplyr
compute_expected_lines <- function(cell_set_meta, cell_line_cols) {
  # Get number of expected cell lines for each cell_set
  result <- cell_set_meta %>%
    dplyr::group_by(across(all_of(cell_line_cols))) %>%
    dplyr::filter(n() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(cell_set) %>%
    dplyr::summarise(
      n_expected_lines = dplyr::n_distinct(across(all_of(cell_line_cols))), # Count unique cell_lines for each cell_set
    ) %>%
    dplyr::ungroup()
  return(result)
}

#' Compute read stats
#'
#' This function calculates summary statistics related to reads and cell line recovery for annotated counts data,
#' combining metadata and performing group-wise computations.
#'
#' @param annotated_counts A data frame containing read data with annotations for barcodes, cell lines, and counts.
#' @param cell_set_meta A data frame containing metadata about cell sets and their expected cell lines.
#' @param group_cols A character vector specifying the grouping columns (default: `c("pcr_plate", "pcr_well")`).
#' @param metric A string indicating the column name of the metric to use for calculations (default: `"n"`).
#' @param cell_line_cols A character vector specifying the column names that define a unique cell line (default: `c("depmap_id", "pool_id")`).
#' @param count_threshold An integer specifying the minimum count threshold for a cell line to be considered recovered (default: `40`).
#'
#' @return A data frame summarizing read statistics for each group, including total reads, expected reads, control barcode reads,
#' recovered cell lines, and their fractions.
#'
#' @import dplyr
compute_read_stats <- function(annotated_counts, cell_set_meta, unknown_counts, group_cols = c("pcr_plate","pcr_well","pert_type"),
                               metric = "n", cell_line_cols = c("depmap_id", "pool_id"), count_threshold = 40) {
  # Compute expected lines from cell_set_meta
  expected_lines <- compute_expected_lines(cell_set_meta, cell_line_cols)

  # Group unknown_counts by group_cols
  unknown_counts <- unknown_counts %>%
      dplyr::left_join(unique(annotated_counts %>% select(pcr_plate, pcr_well, pert_type, pert_plate)),
                       by = c("pcr_plate","pcr_well")) %>%
      dplyr::group_by(across(all_of(group_cols))) %>%
      dplyr::summarise(
      n = sum(.data[[metric]], na.rm = TRUE),
      expected_read = FALSE
      )

  plate_well <- annotated_counts %>%
    dplyr::left_join(expected_lines, by = "cell_set") %>% # Add n_expected_lines from lookup
    # Append unknown_counts, filling NA for any column not present in unknown_counts
    dplyr::bind_rows(unknown_counts) %>%
    dplyr::group_by(across(all_of(group_cols))) %>%
    dplyr::summarise(
      # Total reads
      n_total_reads = sum(.data[[metric]], na.rm = TRUE),
      # Reads mapping to cell lines
      n_expected_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      # Reads mapping to control barcodes
      n_cb_reads = sum(.data[[metric]][cb_name != ""], na.rm = TRUE),
      # Median reads for control barcodes
      median_cb_reads = median(.data[[metric]][cb_name != ""], na.rm = TRUE),
      # Fraction of reads mapping to cell lines
      fraction_expected_reads = n_expected_reads / n_total_reads,
      # Number of cell lines with coverage above 40 reads
      n_lines_recovered = sum(.data[[metric]] >= count_threshold & (is.na(cb_name) | cb_name == ""), na.rm = TRUE),
      # Number of expected lines based on metadata
      n_expected_lines = max(n_expected_lines, na.rm = TRUE), # Bring forward from join
      # Fraction of cell lines with coverage above count threshold
      fraction_cl_recovered = n_lines_recovered / max(n_expected_lines, na.rm = TRUE),
      # Ratio of control barcode reads to cell line reads
      cb_cl_ratio_well = n_cb_reads / n_expected_reads,
    ) %>%
    dplyr::ungroup()

  plate_pert_type <- plate_well %>%
    dplyr::group_by(pcr_plate, pert_type) %>%
    dplyr::summarise(
      cb_cl_ratio_plate = median(cb_cl_ratio_well, na.rm = TRUE),
    ) %>%
    dplyr::ungroup()

    # Combine plate_well and plate_pert_type
  result <- plate_well %>%
    dplyr::left_join(plate_pert_type, by = c("pcr_plate", "pert_type"))

  return(result)
}

#' Calculate CB metrics
#'
#' This function calculates control metrics.
#'
#' @param normalized_counts A data frame containing normalized read counts and associated metrics.
#' @param group_cols A character vector specifying the grouping columns (default: `c("pcr_plate", "pcr_well")`).
#' @param cb_meta A data frame containing control barcode metadata.
#'
#' @return A data frame containing unique combinations of `group_cols` and the calculated metrics.
#'
#' @import dplyr
#'
calculate_cb_metrics <- function(normalized_counts,cb_meta, group_cols = c("pcr_plate", "pcr_well"), pseudocount = 20) {
  valid_profiles= normalized_counts %>% dplyr::filter(!pert_type %in% c(NA, "empty", "", "CB_only"), n != 0,
                                                      cb_ladder %in% unique(cb_meta$cb_ladder),
                                                      cb_name %in% unique(cb_meta$cb_name)) %>%
  dplyr::group_by(across(all_of(group_cols))) %>% dplyr::filter(dplyr::n() > 4) %>% dplyr::ungroup()
  fit_stats= valid_profiles %>%
  dplyr::group_by(across(all_of(group_cols))) %>%
  dplyr::mutate(log2_normalized_n= log2(n+pseudocount) + cb_intercept,
                  cb_mae= median(abs(cb_log2_dose- log2_normalized_n)),
                  mean_y= mean(cb_log2_dose),
                  residual2= (cb_log2_dose- log2_normalized_n)^2,
                  squares2= (cb_log2_dose- mean_y)^2,
                  cb_r2= 1- sum(residual2)/sum(squares2),
                  cb_spearman= cor(cb_log2_dose, log2(n+pseudocount), method= 'spearman', use= 'pairwise.complete.obs')) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(across(all_of(c(group_cols, 'cb_mae', 'cb_r2', 'cb_spearman', 'cb_intercept'))))
  return(fit_stats)
}

# TABLE GENERATION FUNCTION ----------
#' Generate ID column QC metrics table
#'
#' This function computes and combines various quality control (QC) metrics, such as read statistics, skew,
#' and control barcode metrics, grouped by a specified list of id_cols. For standard sequencing runs, this should always
#' be "pcr_plate" and "pcr_well".
#'
#' @param annotated_counts A data frame containing annotated read counts with metadata.
#' @param normalized_counts A data frame containing normalized read counts and associated metrics/metadata.
#' @param cell_set_meta A data frame containing metadata about cell sets and expected cell lines.
#' @param id_cols_list A character vector specifying the columns to group by for QC metric computations.
#' @param cell_line_cols A character vector specifying the column names that define a unique cell line.
#' @param count_threshold An integer specifying the minimum count threshold for a cell line to be considered recovered (default: `40`).
#'
#' @return A data frame (`id_cols_table`) that combines QC metrics, including read statistics, skew, and control barcode metrics, grouped by the specified id_cols.
#'
#' @import dplyr
generate_id_cols_table <- function(annotated_counts, normalized_counts, unknown_counts, cell_set_meta, cb_meta, id_cols_list, cell_line_cols,
                                   count_threshold= 40, pseudocount= 20) {
  paste0("Computing QC metrics grouping by ", paste0(id_cols_list, collapse = ","), ".....")

  read_stats_grouping_cols <- c(id_cols_list, "pert_type", "pert_plate")

  read_stats <- compute_read_stats(annotated_counts = annotated_counts, unknown_counts = unknown_counts, group_cols = read_stats_grouping_cols,
                                   cell_set_meta= cell_set_meta, metric = "n", cell_line_cols = cell_line_cols,
  count_threshold = count_threshold)

  skew <- compute_skew(annotated_counts, group_cols = id_cols_list, metric = "n")

  cb_metrics <- calculate_cb_metrics(normalized_counts, cb_meta, group_cols = id_cols_list, pseudocount = pseudocount)

  id_cols_table <- read_stats %>%
    dplyr::left_join(skew, by = id_cols_list) %>%
    dplyr::left_join(cb_metrics, by = id_cols_list)

  return(id_cols_table)
}

# CELL LINE BY PLATE (PCR_PLATE, CELL LINE) ----------

#' Compute error rate
#'
#' This function calculates the error rate using receiver operating characteristic (ROC) curve data.
#' It quantifies the overlap or misclassification between positive controls and negative controls
#' within grouped data, providing a measure of assay quality. A high error rate indicates large overlap.
#'
#' @param df A data frame containing data for analysis, including metrics and annotations for negative and positive controls.
#' @param metric A string specifying the column name of the metric to use for computing error rate (default: `"log2_normalized_n"`).
#' @param group_cols A character vector specifying the grouping columns (default: `c("depmap_id", "pcr_plate")`).
#' @param negcon A string specifying the `pert_type` value that identifies negative control samples (default: `"ctl_vehicle"`).
#'   These represent baseline or no-treatment conditions and are expected to have lower values.
#' @param poscon A string specifying the `pert_type` value that identifies positive control samples (default: `"trt_poscon"`).
#'   These represent treated conditions and are expected to have higher values.
#'
#' @return A data frame with one row per group and a column `error_rate` containing the computed error rate.
#'
#' @details
#' - **How the Error Rate is Calculated**:
#'   - The function filters the data to include only rows where `pert_type` matches `negcon` or `poscon`.
#'   - For each group defined by `group_cols`, the function computes an ROC curve using the specified `metric`.
#'   - The `metric` values are treated as scores, and the `pert_type` column is used to assign class labels:
#'     - `negcon` samples are treated as the negative class (class 0).
#'     - `poscon` samples are treated as the positive class (class 1).
#'   - The error rate is calculated as the minimum misclassification rate, which is derived from the ROC curve:
#'     - \(\text{Error Rate} = \min\left(\text{True Positive Rate} + 1 - \text{True Negative Rate}\right) / 2\).
#'
#' @import dplyr
#' @import PRROC
compute_error_rate <- function(df, metric = 'log2_normalized_n', group_cols = c("depmap_id", "pcr_plate"),
                               negcon = "ctl_vehicle", poscon = "trt_poscon") {
  paste0("Computing error rate using ", negcon, " and ", poscon, ".....")
  paste0("Grouping by ", paste0(group_cols, collapse = ","), ".....")
  result <- df %>%
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
  return(result)
}

#' Compute control median and MAD
#'
#' This function calculates the median and median absolute deviation (MAD) for negative controls
#' and positive controls in both raw and normalized data. It also computes false sensitivity probabilities
#' based on the MAD values for the negative controls.
#'
#' @param df A data frame containing normalized and raw data along with annotations for control types.
#' @param group_cols A character vector specifying the columns to group by, along with `pert_type` (default: `c("depmap_id", "pcr_plate")`).
#' @param negcon A string specifying the `pert_type` value for negative controls (default: `"ctl_vehicle"`).
#' @param poscon A string specifying the `pert_type` value for positive controls (default: `"trt_poscon"`).
#'
#' @return A data frame containing medians and MADs for raw and normalized data for both control types.
#' Additional columns include false sensitivity probabilities at thresholds -1 (`false_sensitivity_probability_50`)
#' and -2 (`false_sensitivity_probability_25`).
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom stats mad pnorm
compute_ctl_medians_and_mad <- function(df, group_cols = c("depmap_id", "pcr_plate"),
                                    negcon = "ctl_vehicle", poscon = "trt_poscon", pseudocount = 20) {
  paste0("Adding control median and MAD values for ", negcon, " and ", poscon, ".....")
  paste0("Computing falses sensitivity probability for ", negcon, ".....")
  # Group and compute medians/MADs
  result <- df %>%
    dplyr::filter(pert_type %in% c(negcon, poscon)) %>%
    dplyr::group_by(across(all_of(c(group_cols, "pert_type")))) %>%
    dplyr::summarise(
      median_normalized = median(log2_normalized_n, na.rm = TRUE),
      n_replicates = n(),
      mad_normalized = mad(log2_normalized_n, na.rm = TRUE),
      median_raw = median(log2(n+pseudocount), na.rm = TRUE),
      mad_raw = mad(log2(n+pseudocount), na.rm = TRUE)
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
  return(result)
}

#' Compute log fold change of positive controls
#'
#' This function calculates the log fold change (LFC) between positive controls (`poscon`) and negative controls (`negcon`)
#' based on median values from normalized and raw data.
#'
#' @param df A data frame containing median values for both normalized and raw data, with separate columns for positive and negative controls.
#' @param negcon A string specifying the prefix for columns corresponding to the negative control (default: `"ctl_vehicle"`).
#' @param poscon A string specifying the prefix for columns corresponding to the positive control (default: `"trt_poscon"`).
#'
#' @return A data frame with additional columns:
#' - `lfc_normalized`: Log fold change for normalized data.
#' - `lfc_raw`: Log fold change for raw data.
#'
#' @import dplyr
compute_control_lfc <- function(df, negcon = "ctl_vehicle", poscon = "trt_poscon", grouping_cols = c("depmap_id", "pcr_plate")) {
  paste0("Computing log fold change for ", negcon, " and ", poscon, ".....")
  result <- df %>%
    dplyr::mutate(
      lfc_poscon_normalized = .data[[paste0("median_normalized_", poscon)]] -
                       .data[[paste0("median_normalized_", negcon)]],
      lfc_poscon_raw = .data[[paste0("median_raw_", poscon)]] -
                .data[[paste0("median_raw_", negcon)]]
    ) %>%
    dplyr::select(all_of(grouping_cols), lfc_poscon_normalized, lfc_poscon_raw)
  return(result)
}

#' Compute cell line fractions
#'
#' This function computes the total reads and fraction of reads contributed by each cell line within groups defined by specified columns.
#'
#' @param df A data frame containing read data for cell lines, including a metric for read counts.
#' @param metric A string specifying the column name of the metric to use for calculations (default: `"n"`).
#' @param grouping_cols A character vector specifying the columns to group by (default: `c("pcr_plate", "depmap_id")`).
#'
#' @return A data frame with the following additional columns for each group:
#' - `total_reads`: The total reads per group.
#' - `fraction_of_reads`: The fraction of total reads contributed by each entry.
#'
#' @import dplyr
compute_cl_fractions <- function(df, metric = "n", grouping_cols = c("pcr_plate", "depmap_id")) {
  paste0("Computing cell line fractions for ", metric, ".....")
  result <- df %>%
    dplyr::group_by(across(all_of(grouping_cols))) %>%
    dplyr::summarise(
      total_reads = sum(.data[[metric]], na.rm = TRUE),  # Total reads per group
      fraction_of_reads = sum(.data[[metric]], na.rm = TRUE) / sum(.data[[metric]], na.rm = TRUE)  # Fraction of reads for each entry
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(all_of(grouping_cols), total_reads, fraction_of_reads)
  return(result)
}

# TABLE GENERATION FUNCTION ----------

#' Generate cell plate table
#'
#' This function generates a comprehensive QC table for cell_lines +pcr plates by computing and merging various QC metrics, including medians, MADs, error rates, log fold changes (LFC), and cell line fractions.
#'
#' @param normalized_counts A data frame containing normalized read counts and associated metadata.
#' @param filtered_counts A data frame containing filtered read counts for the same dataset.
#' @param cell_line_cols A string of comma-separated column names that define a unique cell line.
#'
#' @return A data frame (`plate_cell_table`) that merges various QC metrics grouped by cell lines and plates, including:
#' - Control medians and MADs for normalized and raw data.
#' - Error rates based on positive and negative control separation.
#' - Log fold changes (LFC) for positive controls.
#' - Fractions of reads contributed by each cell line.
#'
#' @import dplyr
generate_cell_plate_table <- function(normalized_counts, filtered_counts, cell_line_cols, pseudocount = 20) {
  cell_line_list <- strsplit(cell_line_cols, ",")[[1]]
  cell_line_plate_grouping <- c(cell_line_list,"pcr_plate", "pert_plate") # Define columns to group by
  paste0("Computing QC metrics grouping by ", paste0(cell_line_plate_grouping, collapse = ","), ".....")

  # Compute control medians and MAD
  medians_and_mad <- compute_ctl_medians_and_mad(
    df = normalized_counts,
    group_cols = cell_line_plate_grouping,
    negcon = args$negcon_type,
    poscon = args$poscon_type,
    pseudocount = pseudocount
  )

  # Compute error rate
  error_rates <- compute_error_rate(
    df = normalized_counts,
    metric = "log2_normalized_n",
    group_cols = cell_line_plate_grouping,
    negcon = args$negcon_type,
    poscon = args$poscon_type
  )

  # Compute poscon LFC
  poscon_lfc <- compute_control_lfc(
    df = medians_and_mad,
    negcon = args$negcon_type,
    poscon = args$poscon_type,
    grouping_cols = cell_line_plate_grouping
  )

  # Compute cell line fractions per plate
  cell_line_fractions <- compute_cl_fractions(
    df = filtered_counts,
    grouping_cols = cell_line_plate_grouping
  )

  # Merge all tables together
  paste0("Merging ", paste0(cell_line_plate_grouping, collapse = ","), " QC tables together.....")
  plate_cell_table <- medians_and_mad %>%
    dplyr::left_join(error_rates, by = cell_line_plate_grouping) %>%
    dplyr::left_join(poscon_lfc, by = cell_line_plate_grouping) %>%
    dplyr::left_join(cell_line_fractions, by = cell_line_plate_grouping)
  return(plate_cell_table)
}

# QC FLAG FUNCTIONS ----------
id_cols_qc_flags <- function(annotated_counts,
                             normalized_counts,
                             unknown_counts,
                             cb_meta,
                             group_cols = c("pcr_plate","pcr_well","pert_type"),
                             metric = "n",
                             pseudocount = 20,
                             contamination_threshold = 0.8,
                             cb_mae_threshold = 1,
                             cb_spearman_threshold = 0.88,
                             cb_cl_ratio_low_negcon = 0.5,
                             cb_cl_ratio_high_negcon = 2,
                             cb_cl_ratio_low_poscon = 0.5,
                             cb_cl_ratio_high_poscon = 2,
                             well_reads_threshold = 400) {

  # Calculate cb metrics from normalized counts
  cb_metrics <- calculate_cb_metrics(normalized_counts, cb_meta, group_cols = c("pcr_plate", "pcr_well", "pert_type", "pert_plate"), pseudocount = pseudocount)

  # Group unknown_counts by group_cols
  unknown_counts_proc <- unknown_counts %>%
      # Add pert_type from annotated_counts
      dplyr::left_join(unique(annotated_counts %>% select(pcr_plate, pcr_well, pert_type)),
                       by = c("pcr_plate","pcr_well")) %>%
      dplyr::group_by(across(all_of(group_cols))) %>%
      dplyr::summarise(
      n = sum(.data[[metric]], na.rm = TRUE),
      expected_read = FALSE
      )

  # Combine annotated_counts (after adding expected_lines) with unknown_counts and cb_metrics
  combined_data <- annotated_counts %>%
    # Add uknown barcode reads
    bind_rows(unknown_counts_proc) %>%
    # Add cb_metrics
    left_join(cb_metrics, by = group_cols)


  # Initialize a container to record flagged wells.
  flagged_all <- tibble()

  # Working dataset: we will progressively filter out flagged wells.
  working <- combined_data

  ### WELL_READS
  ## Flag wells based on the median count of control barcodes in each well.

  # Calculate the median control barcode reads for each well.
  well_reads <- working %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_cb_reads        = sum(.data[[metric]][cb_name != ""], na.rm = TRUE),
      median_cb_reads   = median(.data[[metric]][cb_name != ""], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Flag wells where median_cb_reads is less than 20 * pseudocount.
    mutate(qc_flag = if_else(median_cb_reads < (well_reads_threshold), "well_reads", NA_character_))

  # Record flagged wells and remove them from the working dataset.
  flagged_well_reads <- well_reads %>%
    filter(!is.na(qc_flag)) %>%
    select(all_of(group_cols), qc_flag)
  flagged_all <- bind_rows(flagged_all, flagged_well_reads)

  working <- working %>% anti_join(flagged_well_reads, by = group_cols)

  ### CONTAMINATION
  ## Flag wells based on the fraction of reads mapping to expected cell lines or control barcodes.

  # Calculate the total and expected reads for each well.
  containation <- working %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_total_reads    = sum(.data[[metric]], na.rm = TRUE),
      n_expected_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(fraction_expected_reads = n_expected_reads / n_total_reads,
           qc_flag = if_else(fraction_expected_reads < contamination_threshold, "contamination", NA_character_))

  # Record flagged wells and remove them from the working dataset.
  flagged_contamination <- containation %>%
    filter(qc_flag == "contamination") %>%
    select(all_of(group_cols), qc_flag)
  flagged_all <- bind_rows(flagged_all, flagged_contamination)

  working <- working %>% anti_join(flagged_contamination, by = group_cols)

  ### CB_LINEARITY
  ## Flag wells based on the linearity of the control barcodes in each well defined by the median absolute error (MAE) and Spearman correlation coefficient.

  # Calculate the MAE and Spearman correlation
  cb_linearity <- working %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      cb_mae = median(.data[["cb_mae"]], na.rm = TRUE),
      cb_spearman = median(.data[["cb_spearman"]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(qc_flag = if_else(cb_mae > cb_mae_threshold | cb_spearman < cb_spearman_threshold, "cb_linearity", NA_character_))

  # Record flagged wells and remove them from the working dataset.
  flagged_cb_linearity <- cb_linearity %>%
    filter(qc_flag == "cb_linearity") %>%
    select(all_of(group_cols), qc_flag)
  flagged_all <- bind_rows(flagged_all, flagged_cb_linearity)
  working <- working %>% anti_join(flagged_cb_linearity, by = group_cols)

  ### CB_CL_RATIO
  ## POSITIVE CONTROL
  ## Flag wells based on the ratio of control barcode reads to cell line reads.
  cb_cl_ratio_poscon = working %>%
    filter(pert_type == "trt_poscon") %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_expected_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      n_cb_reads = sum(.data[[metric]][cb_name != ""], na.rm = TRUE),
      cb_cl_ratio_well = n_cb_reads / n_expected_reads,
      .groups = "drop"
    ) %>%
    group_by(pcr_plate, pert_type) %>%
    mutate(cb_cl_ratio_plate = median(cb_cl_ratio_well, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(qc_flag = if_else(cb_cl_ratio_well < cb_cl_ratio_low_poscon * cb_cl_ratio_plate | cb_cl_ratio_well > cb_cl_ratio_high_poscon * cb_cl_ratio_plate, "cb_cl_ratio", NA_character_))

  # Record flagged wells and remove them from the working dataset.
  flagged_cb_cl_ratio_poscon <- cb_cl_ratio_poscon %>%
    filter(qc_flag == "cb_cl_ratio") %>%
    select(all_of(group_cols), qc_flag)
  flagged_all <- bind_rows(flagged_all, flagged_cb_cl_ratio_poscon)
  working <- working %>% anti_join(flagged_cb_cl_ratio_poscon, by = group_cols)

  ## NEGATIVE CONTROL
  ## Flag wells based on the ratio of control barcode reads to cell line reads.
  cb_cl_ratio_negcon = working %>%
    filter(pert_type == "ctl_vehicle") %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_expected_reads = sum(.data[[metric]][expected_read], na.rm = TRUE),
      n_cb_reads = sum(.data[[metric]][cb_name != ""], na.rm = TRUE),
      cb_cl_ratio_well = n_cb_reads / n_expected_reads,
      .groups = "drop"
    ) %>%
    group_by(pcr_plate, pert_type) %>%
    mutate(cb_cl_ratio_plate = median(cb_cl_ratio_well, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(qc_flag = if_else(cb_cl_ratio_well < cb_cl_ratio_low_negcon * cb_cl_ratio_plate | cb_cl_ratio_well > cb_cl_ratio_high_negcon * cb_cl_ratio_plate, "cb_cl_ratio", NA_character_))

  # Record flagged wells and remove them from the working dataset.
  flagged_cb_cl_ratio_negcon <- cb_cl_ratio_negcon %>%
    filter(qc_flag == "cb_cl_ratio") %>%
    select(all_of(group_cols), qc_flag)
  flagged_all <- bind_rows(flagged_all, flagged_cb_cl_ratio_negcon)


  ### RETURN RESULTS
  # Filter normalized_counts for only the wells that were not flagged
  normalized_filtered <- flagged_all %>%
    select(pcr_plate, pcr_well, pert_type) %>%
    unique() %>%
    dplyr::left_join(normalized_counts, by = group_cols)

  # Return a list containing both the filtered normalized_counts and a record of all flagged wells.
  list(result = normalized_filtered, well_flags = flagged_all)
}

## POOL WELL QC FLAGS
pool_well_qc_flags <- function(normalized_counts, pool_well_qc) {
  ### POOL_WELL_OUTLIERS
  ## Flag pool/well combinations based on the fraction of cell lines in a pool + well that are some distance from the pool + well median.
  ## This is in control wells only
  flagged_all <- pool_well_qc %>%
    filter(!is.na(qc_flag))

  ### RETURN RESULTS
  # Filter normalized_counts for only the wells that were not flagged
  normalized_filtered <- normalized_counts %>%
    dplyr::left_join(flagged_all %>%
      select(pcr_plate, pcr_well, pert_type, pool_id),
                     by = c("pool_id", "pcr_plate", "pcr_well", "pert_type"))

  # Return a list containing both the filtered normalized_counts and a record of all flagged wells.
  list(result = normalized_filtered, pool_well_flags = flagged_all)
}

generate_pool_well_qc_table <- function(normalized_counts, pool_well_delta_threshold = 5, pool_well_fraction_threshold = 0.4) {
  ### POOL_WELL_OUTLIERS
  ## Flag pool/well combinations based on the fraction of cell lines in a pool + well that are some distance from the pool + well median.
  ## This is in control wells only
  pool_well_outliers <- normalized_counts %>%
    # Consider only cell line barcodes
    filter(is.na(cb_name)) %>%
    # Get the median value of the pool in each well
    group_by(pcr_plate, pcr_well, pert_type, pool_id) %>%
    mutate(
      pool_well_median = median(log2_normalized_n, na.rm = TRUE),
      n_cell_lines = n_distinct(paste(lua, depmap_id, cell_set, sep = "_"))
    ) %>%
    mutate(
      delta_from_pool_well_median = abs(log2_normalized_n - pool_well_median)
    ) %>%
    group_by(pcr_plate, pcr_well, pert_type, pool_id) %>%
    reframe(
      n_outliers = sum(delta_from_pool_well_median > pool_well_delta_threshold, na.rm = TRUE),
        fraction_outliers = n_outliers / n_cell_lines
    ) %>%
    ungroup() %>%
    mutate(qc_flag = if_else(fraction_outliers > pool_well_fraction_threshold, "pool_well_outliers", NA_character_))
}

