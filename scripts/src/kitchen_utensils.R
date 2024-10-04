# Kitchen Utensils - 
# This file contains functions for the pipeline.
# The functions are sorted alphabetically

#' process_in_chunks
#' 
#' This function runs some action over chunks of a large file. At the end, returns a list of all the chunks
#' 
#' @param large_file_path description
#' @param chunk_size description
#' @param action A function passed to act on each chunk
#' @param ... Additional parameters to be passed into the action parameter
process_in_chunks= function(large_file_path, chunk_size= 10^6, action, ...) {
  header_col_names= data.table::fread(large_file_path, header= TRUE, sep= ',', nrow= 0) %>% colnames()
  chunk_idx= 1 # Counter to keep track of chunks in a loop
  current_chunk_size= chunk_size # Variable for loop exit condition
  chunk_collector= list() # List to collect processed chunks
  
  # For each chunk, call an action
  while(current_chunk_size == chunk_size) {
    current_chunk= data.table::fread(large_file_path, header= FALSE, sep= ',',
                                     col.names= header_col_names,
                                     nrow= chunk_size, skip= chunk_size * (chunk_idx - 1) + 1)
    
    current_chunk_size= nrow(current_chunk) # set current chunk size to stop loop
    print(paste('Working on chunk', chunk_idx, 'with', current_chunk_size, 'rows.', sep= ' '))
    
    chunk_collector[[chunk_idx]]= do.call(action, list(current_chunk, ...))
    chunk_idx= chunk_idx + 1
  }
  
  # Return a list of all the chunks
  return(chunk_collector)
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

#' validate_unique_samples
#' 
#' This function checks that a list of columns uniquely identifies every row of a dataframe.
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