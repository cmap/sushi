#' validate_columns_exist
#' 
#' This function checks that a list of columns are present in a dataframe.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_columns_exist= function(selected_columns, df) {
  # Check that all of selected_columns are in df
  if(any(!selected_columns %in% colnames(df))) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
  
#' validate_columns_entries
#' 
#' This function checks that for a list of columns, all entries are filled in.
#' It checks all column entries against a list of potential empty values.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @param empty_values Optional vector of values that equate to empty. Defaults to NA, "NA", "", and " ".
#' @return Boolean
validate_columns_entries= function(selected_columns, df, empty_values= c(NA, 'NA', '', ' ')) {
  # Check for rows in selected_columns that equate to predefined empty values.
  missing_rows= df %>% dplyr::filter(if_any(all_of(selected_columns), ~ . %in% empty_values))
  if(nrow(missing_rows) > 0) {
    print('The following rows in the sample meta are not filled out for the sequencing index columns.')
    print(missing_rows) # show the empty rows
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' validate_unique_samples
#' 
#' This function checks that a list of columns uniquely identifies all entries of a dataframe.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_unique_samples= function(selected_columns, df) {
  unique_column_values= df %>% dplyr::distinct(pick(all_of(selected_columns)))
  if(nrow(unique_column_values) != nrow(df)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' collate_fastq_reads
#' 
#' This function takes in the fastq reads (uncollapsed_raw_counts) and
#' filters for reads coming from flowcells specificed in the sample meta.
#' The function then sums up the reads across specified sequencing index columns.
#' 
#' @param uncollapsed_raw_counts Data frame of reads from all the fastq files with the following columns - \cr
#'                    "flowcell_name", "flowcell_lane", "index_1", "index_2", and "forward_read_cl_barcode", "n"
#' @param sample_meta Sample metadata generate for the project which may contain the following columns - 
#'                    "flowcell_names", "flowcell_lanes", "index_1", "index_2". The sample meta must contain
#'                    "flowcell_names" and "flowcell_lanes" for filtering.
#' @param sequencing_index_cols Sequencing columns from the sample meta that the counts should be collapsed on. \cr
#'                              This defaults onto the following columns: "index_1", "index_2"
#' @returns Returns a dataframe with columns specified by the sequencing_index_cols, "forward_read_cl_barcode", and "n".
#' @import tidyverse
collate_fastq_reads= function(uncollapsed_raw_counts, sample_meta, sequencing_index_cols= c('index_1', 'index_2')) {
  require(tidyverse)
  
  # Validation: Check that flowcell_names and flowcell_lanes exist in the sample meta ----
  if(!validate_columns_exist(c('flowcell_names', 'flowcell_lanes'), sample_meta)) {
    stop('flowcell_names and/or flowcell_lanes is NOT present in the sample meta.')
  }
  
  # Validation: Check that sequencing_index_cols exist in the sample meta ----
  if(!validate_columns_exist(sequencing_index_cols, sample_meta)) {
    stop('One or more sequencing_index_cols is NOT present in the sample meta.')
  }
  
  # Validation: Check that sequencing_index_cols in the sample meta are filled out ----
  # Check for rows in sequencing_index_cols that equate to empty - NA, "NA", "", " "
  # Error out of the sequencing_index_cols are not filled out in the sample meta.
  if(!validate_columns_entries(sequencing_index_cols, sample_meta)) {
    stop('One or more sequencing_index_cols in the sample meta is not filled out.')
  }
  
  # Validation: Check that sequencing_index_cols uniquely identify rows of sample meta ----
  if(!validate_unique_samples(sequencing_index_cols, sample_meta)) {
    print('There may be multiple entries in the sample meta that have the same combination of sequencing index columns.')
    stop('The specified sequencing index columns do NOT uniquely identify every PCR well.')
  }
  
  # Determine which flowcell names + lanes are expected ----
  # "flowcell_names" and "flowcell_lanes" are strings that can contain more than one item.
  # Columns can be parsed by splitting on the chars , ; :
  # If there are multiple lane names and lane numbers, this uses the Cartesian product!
  # Note: fread and read.csv keeps commas, read_csv DROPS commas
  meta_flowcells= sample_meta %>% dplyr::distinct(flowcell_names, flowcell_lanes) %>%
    dplyr::mutate(flowcell_name= base::strsplit(flowcell_names, split='[,;:]', fixed=F),
                  flowcell_lane= base::strsplit(flowcell_lanes, split='[,;:]', fixed=F)) %>% 
    tidyr::unnest(cols= flowcell_name) %>% tidyr::unnest(cols= flowcell_lane) %>%
    dplyr::mutate(flowcell_lane= as.numeric(flowcell_lane))
  
  # Print out expected flowcells from the sample meta.
  print(paste0('Identified ', nrow(meta_flowcells), ' unique flowcell + lane combos in the sample meta ...'))
  print(meta_flowcells)
  
  # Print warning if there are multiple flowcell names with multiple flowcell lanes.
  multi_name_and_lanes= meta_flowcells %>% dplyr::filter(grepl(',:;', flowcell_names) & grepl(',:;', flowcell_names))
  if(nrow(multi_name_and_lanes) > 0) {
    print('WARNING: Detected sample(s) sequenced over multiple flowcells and flowcell lanes.')
    print('The function assumes that the same lanes were used for both flowcells.')
  }
  
  # Validation: Check that all expected flowcell name + lanes are detected ----
  # Check that all expected flowcell name + lanes are present in uncollapsed raw counts.
  missing_flowcells= meta_flowcells %>%
    dplyr::anti_join(uncollapsed_raw_counts %>% dplyr::distinct(flowcell_name, flowcell_lane),
                     by= c('flowcell_name', 'flowcell_lane'))
  if(nrow(missing_flowcells)!= 0) {
    print('The following flowcells/lanes specified in the sample meta were not detected in the fastq reads.')
    print(missing_flowcells)
    print('Check that the sample meta is correct or that all fastq files are in the correct directory.')
    stop('One or more flowcell specified in the sample meta was not detected.')
  }
  
  # Create raw counts by summing over appropriate sequencing_index_cols ----
  # Use an inner join to collect reads with valid flowcell name/lane combinations, 
  # then sum reads across sequencing columns
  raw_counts= uncollapsed_raw_counts %>% 
    dplyr::inner_join(meta_flowcells, by= c('flowcell_name', 'flowcell_lane'), relationship= 'many-to-one') %>%
    dplyr::group_by(pick(all_of(c(sequencing_index_cols, 'forward_read_cl_barcode')))) %>% 
    dplyr::summarize(n= sum(n)) %>% dplyr::ungroup()
  
  return(raw_counts)
}
