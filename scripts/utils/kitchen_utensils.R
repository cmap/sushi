# Kitchen Utensils -
# This file contains functions for the pipeline.
# The functions are sorted alphabetically

#' process_in_chunks
#'
#' This function runs some action over chunks of a large file. At the end, returns a list of all the chunks
#'
#' @param large_file_path Path to a large csv file. This file may be too large to read into R.
#' @param chunk_size The number of rows in a chunk.
#' @param action A function to perform over a chunk.
#' @param ... Additional parameters to be passed into the action parameter
process_in_chunks= function(large_file_path, chunk_size= 10^6, action, ...) {
  # Read in the column names. These names will be passed onto each chunk.
  # When reading a file in chunks, the column names in the first line are not always passed.
  # Use data.table to read in just the headers with nrow= 0.
  header_col_names= read_data_table(large_file_path, header= TRUE, sep= ',', nrow= 0) %>% colnames()
  chunk_idx= 1 # Counter to keep track of chunks in a loop
  current_chunk_size= chunk_size # Variable for loop exit condition
  chunk_collector= list() # List to collect processed chunks

  # For each chunk, call an action
  while(current_chunk_size == chunk_size) {
    # Read in a chunk of the large file and set the column names.
    # nrow - the number of rows to read in
    # skip - the number of rows to skip before starting to read in.
    current_chunk= read_data_table(large_file_path, header= FALSE, sep= ',',
                                     col.names= header_col_names,
                                     nrow= chunk_size, skip= chunk_size * (chunk_idx - 1) + 1)

    current_chunk_size= nrow(current_chunk) # set current chunk size to stop loop
    print(paste('Working on chunk', chunk_idx, 'with', current_chunk_size, 'rows.', sep= ' '))

    # Call the action over the chunk
    chunk_collector[[chunk_idx]] <- do.call(
    action,
    c(list(uncollapsed_raw_counts = current_chunk), list(...))
    )

    chunk_idx= chunk_idx + 1
  }

  # Return a list of all the chunks
  return(chunk_collector)
}

#' Read a CSV file with enforced data types from a master schema using data.table
#'
#' This function reads a CSV file into a data.table, applying column types
#' specified in a central schema.yaml file. It is designed for high performance
#' by leveraging read_data_table. It uses the 'here' package to reliably
#' locate the schema file relative to the project root.
#'
#' @param csv_path The path to the CSV file to be read.
#' @param schema_path The path to the schema YAML file. Defaults to
#'   a path constructed from the project root.
#'
#' @return A data.table with data types enforced according to the schema.
#'
#' @import data.table
#' @import yaml
#' @import here
#'
#' @examples
#' \dontrun{
#' df <- read_csv_with_schema("path/to/your/data.csv")
#' }
read_data_table <- function(csv_path, schema_path = NULL) {

  # Ensure required packages are available
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed.")
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required but not installed.")
  }

  # --- Define the mapping from abstract types to data.table class names ---
  # Note: For performance, datetimes are read as character and converted after.
  type_mapping <- c(
    "integer" = "integer",
    "float" = "numeric",
    "string" = "character",
    "boolean" = "logical",
    "datetime" = "character" # Read as character first for speed
  )

  # --- Determine the schema path ---
  if (is.null(schema_path)) {
    # Automatically find the schema file from the project root
    schema_path <- here::here("sushilib", "io", "schema.yaml")
  }

  if (!file.exists(schema_path)) {
    stop(paste("Schema file not found at the expected path:", schema_path))
  }

  # --- Read the master schema file ---
  schema <- yaml::read_yaml(schema_path)

  # --- Create the colClasses vector for fread ---
  # This creates a named character vector, e.g., c(col_a = "integer", col_b = "character")
  col_classes <- vapply(schema, function(dtype) type_mapping[[dtype]], FUN.VALUE = character(1))

  cat("Applying data.table Schema:\n")
  print(col_classes)

  # --- Read the CSV using fread ---
  dt <- read_data_table(csv_path, colClasses = col_classes)

  # --- Post-processing for datetime columns ---
  # Find which columns were originally specified as 'datetime'
  datetime_cols <- names(schema)[which(unlist(schema) == "datetime")]

  if (length(datetime_cols) > 0) {
    cat("\nConverting datetime columns:", paste(datetime_cols, collapse = ", "), "\n")
    # Use data.table's efficient `:=` to convert columns by reference
    for (col in datetime_cols) {
      dt[, (col) := as.POSIXct(get(col))]
    }
  }

  return(dt)
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

#' check_file_exists
#'
#' This function checks that a given file exists, and if not, stops execution.
#' This is meant to be used at the end of a module to validate that outputs have been generated.
#'
#' @param file_path Path to a file to check
check_file_exists <- function(file_path) {
  if (file.exists(file_path)) {
    message(paste("File ", file_path, " was successfully generated."))
  } else {
    stop(paste("Error: File ", file_path, " does not exist, please check the output of this module."))
  }
}

#' Filter out control barcodes
#'
#' This function filters out control barcodes from a dataframe.
#' tryCatch is used in the event that the `cb_name` column does not exist, it then returns the original dataframe.
#'
#' @param df A dataframe to filter.
#' @return A filtered dataframe.
filter_control_barcodes <- function(df) {
  df %>% dplyr::filter(tryCatch(is.na(cb_name), error = function(e) TRUE))
}

#' Append a print statement to a file to track critical console output
#'
#' This function, when applied to a given print statement, will append the statement to a file in order to
#' track critical console outputs for human review. It will also write the statement to the console as usual.
#'
#' @param print_statement A string to print and append to a file
append_critical_output <- function(statement, out) {
  # If the out/logs directory does not exist, create it
    if (!dir.exists(paste0(out, "/logs"))) {
        dir.create(paste0(out, "/logs"), recursive = TRUE)
    }
  # Convert data frames or lists into a string if necessary
  if (!is.character(statement)) {
    statement <- capture.output(print(statement))
  }
  # Print to console
  cat(statement, sep = "\n")
  # Append to file
  cat(statement, file = paste0(out, "/logs/critical_output.txt"), append = TRUE, sep = "\n")
}

#' Delete existing files when re-running a module
#'
#' This function deletes existing files in a given output directory that match a given pattern.
#' This is useful when re-running a module to ensure that old files do not interfere with
#' the new output.
#' @param out_dir The output directory to check for existing files.
#' @param pattern A pattern to match files against. This can be a regular expression.
delete_existing_files <- function(out_dir, pattern) {
  files_to_delete <- list.files(path = out_dir, pattern = pattern, full.names = TRUE)
  if (length(files_to_delete) > 0) {
    file.remove(files_to_delete)
    message(paste("Deleted", length(files_to_delete), "existing files in", out_dir))
    message("Deleted files:")
    print(files_to_delete)
  } else {
    message(paste("No existing files to delete in", out_dir))
  }
}
