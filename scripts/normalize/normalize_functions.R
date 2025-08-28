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
normalize <- function(X, id_cols, CB_meta, pseudocount) {
  # Required functions
  require(magrittr)
  require(tidyverse)

  # Create log2_n with pseudocount ----
  X %<>% dplyr::mutate(log2_n = log2(n + pseudocount))

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
  # Drop wells with invalid pert_type, wells without control barcodes, cell line entries or other CBs,
  # cbs with zero reads, and profiles with fewer than 4 CBs.
  valid_profiles= X %>% dplyr::filter(!pert_type %in% c(NA, "empty", "", "CB_only"), n != 0,
                                      cb_ladder %in% unique(CB_meta$cb_ladder),
                                      cb_name %in% unique(CB_meta$cb_name)) %>%
    dplyr::group_by(dplyr::pick(tidyselect::all_of(id_cols))) %>%
    dplyr::filter(dplyr::n() > 4) %>% dplyr::ungroup()

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
    dplyr::select(-log2_n)

  return(normalized)
}
