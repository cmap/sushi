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

#' join_sample_meta
#' 
#' Joins a given data frame with the sample meta.
#' 
#' @param df
#' @param sample_meta Dataframe of the sample meta used in the run
#' @param key_cols Vector of column names used as identifiers in the sample meta.
#' @returns Data frame with additional columns from the sample meta
join_sample_meta= function(df, sample_meta, key_cols) {
  # Validation: Check that key_cols are present in df ----
  if(validate_columns_exist(key_cols, df) == FALSE) {
    print(key_cols)
    stop('Not all key_cols (printed above) are present in the provided data frame.')
  }
  
  # Validation: Check that key_cols are present in the sample meta ----
  if(validate_columns_exist(key_cols, sample_meta) == FALSE) {
    print(key_cols)
    stop('Not all key_cols (printed above) are present in the sample meta.')
  }
  
  # Collapse the sample meta using key_cols and join onto the input df ----
  collapsed_metadata= sample_meta %>% dplyr::group_by(pick(all_of(key_cols))) %>% 
    dplyr::summarise(across(everything(), function(x) paste(sort(unique(x)), collapse= ' | '))) %>% dplyr::ungroup()

  expanded_df= df %>% dplyr::left_join(collapsed_metadata, by= key_cols, relationship='many-to-one')
  
  # Validation: Check for duplicate columns ----
  duplicate_columns= setdiff(c(colnames(df), colnames(sample_meta)), colnames(expanded_df))
  if(length(duplicate_columns > 0)) {
    print("WARNING: The following column(s) appear in the dataframe and the sample meta, but not in key_cols.")
    print(duplicate_columns)
    print('The columns(s) thus appear twice in the output dataframe.')
  }
  
  return(collapsed_counts)
}