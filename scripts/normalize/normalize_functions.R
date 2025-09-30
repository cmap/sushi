#'  normalize
#'
#'  takes a filtered dataframe of raw read counts and normalizes
#'  counts using control barcodes
#'
#' @param X A dataframe of annotated readcounts that must include the following columns:
#'           log2_n or n: raw readcounts or log2(n) of read counts. Computes log2_n if not present
#'           cb_ladder: column indicating which CB ladder was used if any.
#'           cb_log2_dose: log2 of dose at which control barcode was spiked in, if applicable
#'           cb_name: contains the name of the control barcode that the read corresponds to, or NA
#' @param id_cols a vector of columns used to identify each PCR well and to be carried forward in the pipeline.
#'                These column names should be present in dataframe X.
#' @param CB_meta CB_meta dataframe with the required columns "cb_ladder" and "cb_name".
#' @param pseudocount A pseudocount to be added to the counts so that logs can be taken
#' @returns Dataframe with counts normalized to control barcodes
#' @import tidyverse
#' @import magrittr
normalize <- function(X, id_cols, CB_meta, req_negcon_reps, negcon_type = "ctl_vehicle",
                      pseudocount = 0, cb_mad_cutoff = 1) {
  # Required functions
  require(magrittr)
  require(tidyverse)

  # Create log2_n with pseudocount ----
  X %<>% dplyr::mutate(log2_n = log2(n + pseudocount))

  # Filter out any duplicate cell lines if pool_id has a particular string ----
  if('pool_id' %in% colnames(X)) {
    X %<>% dplyr::filter(!grepl("_+_", pool_id))
  }

  # Validation: Check that id_cols are present in the dataframe ----
  if(!validate_columns_exist(id_cols, X)) {
    stop('One or more id_cols (printed above) is NOT present in the supplied dataframe.')
  }

  # Filter cb_meta by dropping any control barcodes without "well_norm"
  # indicated under "cb_type"
  if ("cb_type" %in% colnames(CB_meta)) {
    dropped_cbs = CB_meta |> dplyr::filter(cb_type != "well_norm")

    if (nrow(dropped_cbs) > 0) {
      print(" The following CBs are excluded from normalization.")
      print(dropped_cbs)
      CB_meta = CB_meta |> dplyr::filter(cb_type == "well_norm")
    }
  }

  # Identify valid profiles and valid control barcodes to determine intercept ----
  # Valid CBs per plate
  message("Calculating MADs for the control barcodes...")
  cb_mad = get_cb_mad(X[!is.na(cb_name) & pert_type == negcon_type],
                      id_cols = id_cols, cb_mad_cutoff = cb_mad_cutoff)

  print("The following CBs are dropped due to high MAD.")
  invalid_cbs = cb_mad |> dplyr::filter(keep_cb == FALSE, num_reps >= req_negcon_reps)
  print(invalid_cbs)

  print("The MAD filter did not apply to these CBs due to low replicate counts.")
  print(cb_mad |> dplyr::filter(keep_cb == FALSE, num_reps < req_negcon_reps))

  # Drop wells with invalid pert_type, wells without control barcodes, cell line entries or other CBs,
  # cbs with zero reads, cbs with high MADs, and profiles with fewer than 4 CBs.
  valid_profiles= X %>% dplyr::filter(!pert_type %in% c(NA, "empty", "", "CB_only"), n != 0,
                                      cb_ladder %in% unique(CB_meta$cb_ladder),
                                      cb_name %in% unique(CB_meta$cb_name)) %>%
    dplyr::group_by(dplyr::pick(tidyselect::all_of(id_cols))) %>%
    dplyr::filter(dplyr::n() > 4) %>% dplyr::ungroup()

  # Drop high MAD CBs
  if (nrow(invalid_cbs) > 0) {
    valid_profiles = valid_profiles |>
      dplyr::anti_join(invalid_cbs, by = c("pcr_plate", "cb_ladder", "cb_name"))
  }

  # Validation: Check which wells/profiles were dropped ----
  distinct_all_profiles= X %>% dplyr::distinct(dplyr::pick(tidyselect::all_of(id_cols)))
  distinct_valid_profiles= valid_profiles %>% dplyr::distinct(dplyr::pick(tidyselect::all_of(id_cols)))
  if(nrow(distinct_all_profiles) != nrow(distinct_valid_profiles)) {
    # Print error if all profiles were dropped
    if(nrow(valid_profiles) == 0) {
      stop('No valid profiles detected for normalization!')
    }

    # Print out the profiles that were dropped
    profiles_dropped_at_norm= distinct_all_profiles %>% dplyr::anti_join(distinct_valid_profiles, by= id_cols)
    print(paste('Number of profiles dropped at normalization:', nrow(profiles_dropped_at_norm),
                'out of', nrow(distinct_valid_profiles)))
    print('Reason for dropping:\n1. trt type is empty, NA, or CB_only. 2. Detected <=4 CBs.')
    print("Dropped profiles are ...")
    print(profiles_dropped_at_norm)
  }

  # Calculate fit intercept for valid profiles using median intercept ----
  fit_intercepts= valid_profiles %>% dplyr::group_by(pick(all_of(c(id_cols, 'cb_log2_dose')))) %>%
    dplyr::summarize(dose_intercept= mean(cb_log2_dose) - mean(log2_n)) %>%
    dplyr::group_by(pick(all_of(id_cols))) %>%
    dplyr::summarize(cb_intercept= median(dose_intercept)) %>% dplyr::ungroup()

  # Normalize entries ----
  normalized= X %>% dplyr::inner_join(fit_intercepts, by=id_cols) %>%
    dplyr::mutate(log2_normalized_n= log2_n + cb_intercept) %>%
    dplyr::left_join(cb_mad, by = c("pcr_plate", "cb_ladder", "cb_name")) %>%
    dplyr::mutate(cb_ladder = ifelse(keep_cb == FALSE, paste0(cb_ladder, " - dropped"), cb_ladder)) |>
    dplyr::select(-log2_n, -keep_cb)

  return(normalized)
}

#' get_cb_mad
#'
#' Calculates MAD for each control barcodes across all negative control wells of a plate
#'
#' @param cb_reads Dataframe of control barcode read counts.
#' @param id_cols Vector of columns that uniquely identify each PCR well.
#' @param cb_mad_cutoff Maximum MAD value for a control barcode.
get_cb_mad = function(cb_reads, id_cols, cb_mad_cutoff = 1) {
  cb_mad = cb_reads |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(id_cols))) |>
    dplyr::mutate(cb_frac = (n + 1) / sum(n)) |>
    dplyr::group_by(pcr_plate, cb_ladder, cb_name) |>
    dplyr::summarise(mad_log2_cb_frac = mad(log2(cb_frac)),
                     num_reps = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(keep_cb = ifelse(mad_log2_cb_frac < cb_mad_cutoff, TRUE, FALSE))

  return(cb_mad)
}

#' add_pseudovalues
#'
#' After normalization without a pseudocount, compute a pseudovalue using the negative controls.
#' Add the pseudovalue to normalized counts.
#'
#' @param norm_counts Dataframe of normalized counts.
#' @param negcon_cols Vector of columns names in norm_counts that describe a negative control condition.
#' @param read_detection_limit Read count detection limit.
#' @param negcon_type String in column pert_type of norm_counts that indicates negative control samples.
#' @returns Dataframe with normalized count adjusted with a pseudovalue.
#' @import tidyverse
add_pseudovalue = function(norm_counts, negcon_cols, read_detection_limit = 10, negcon_type = "ctl_vehicle") {
  # Pseudovalue is calculated over negative control groups.
  # Check that all of those columns are present.
  missing_cols = setdiff(negcon_cols, colnames(norm_counts))
  if (length(missing_cols) > 0) {
    stop("normalize - add_pseudovalue: The following column(s) are missng from normalized_counts: ",
         paste(missing_cols, collapse = ","))
  }

  # Calculate the pseudovalue over each control group.
  pv_per_group = norm_counts |> dplyr::filter(pert_type == negcon_type) |>
    dplyr::group_by(dplyr::across(tidyselect::all_of(negcon_cols))) |>
    dplyr::summarise(log2_pseudovalue = median(log2(read_detection_limit) + cb_intercept), .groups = "drop")

  # Join pseudovalues to the main df and add them to all rows/samples.
  norm_counts = norm_counts |>
    dplyr::left_join(pv_per_group, by = negcon_cols) |>
    dplyr::mutate(log2_normalized_n = log2(2^log2_normalized_n + 2^log2_pseudovalue))

  return(norm_counts)
}