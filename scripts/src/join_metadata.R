#' validate_columns_exist
#' 
#' This function checks that a list of columns are present in a dataframe.
#' Columns that were not found in the dataframe are printed out.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_columns_exist= function(selected_cols, df) {
  # Check that all of selected_columns are in df - base::setdiff(A, B) = A[!A %in% B].
  unmatched_cols= base::setdiff(selected_cols, colnames(df))
  
  if(length(unmatched_cols) > 0) {
    print('The following columns are missing: ')
    print(unmatched_cols)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Join metadata
#' 
#' Joins a given data frame with the sample meta.
#' 
#' @param input_df Input dataframe that should contain the columns specified in the "key_cols" parameter and "cell_set".
#' @param metadata Dataframe of the sample meta used in the run.
#' @param key_cols Vector of column names used as identifiers in the sample meta.
#' @returns Data frame with additional columns from the sample meta.
join_metadata= function(input_df, metadata, key_cols) {
  # Validation: Check that key_cols are present in df ----
  if(validate_columns_exist(key_cols, input_df) == FALSE) {
    stop('Not all key_cols (printed above) are present in the provided dataframe.')
  }
  
  # Validation: Check that key_cols are present in the sample meta ----
  if(validate_columns_exist(key_cols, metadata) == FALSE) {
    stop('Not all key_cols (printed above) are present in the provided metadata.')
  }
  
  # Collapse the sample meta using key_cols and join onto the input df ----
  # Collapse unique values into a single row and then filter out columns with the separator.
  # Columns with only one unique value in a group are selected.
  collapsed_metadata= metadata %>% dplyr::group_by(pick(all_of(key_cols))) %>% 
    dplyr::summarise(across(everything(), function(x) paste(sort(unique(x)), collapse= ':::'))) %>% dplyr::ungroup() %>%
    dplyr::select(all_of(key_cols), where(function(x) base::any(!grepl(':::', x))))
  
  # Join using the key_cols, drop any columns that were duplicated.
  output_df= dplyr::left_join(input_df, collapsed_metadata, by= key_cols, 
                              suffix= c('', '.y'), relationship='many-to-one') %>%
    dplyr::select(-tidyselect::ends_with('.y'))
  
  # Validation: Check that merge did not explode ----
  print(paste0(' Input df rows: ', nrow(input_df)))
  print(paste0('Output df rows: ', nrow(output_df)))
  if(nrow(input_df) < nrow(output_df)) {
    stop('Metadata join is producing more rows than expected!')
  } else if(nrow(input_df) > nrow(output_df)) {
    stop('Metadata join is dropping some rows!')
  } else {}
  
  # Print out the sample meta columns that were added to the dataframe ----
  added_cols= base::setdiff(colnames(output_df), colnames(input_df))
  if(length(added_cols > 0)) {
    print(paste0('The following ', length(added_cols), ' column(s) were added:'))
    print(added_cols)
    print(paste0('The following ', length(metadata) - length(added_cols) - length(key_cols),
                 ' column(s) from the metadata were not added. They may already exist in the dataframe.'))
    print(base::setdiff(colnames(metadata), c(added_cols, key_cols)))
  } else {
    print('No additional columns from the metadata were added.')
  }
  
  return(output_df)
}
