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

#' validate_detected_flowcells
#' 
#' This function checks that all the expected flowcells are present in a table of detected flowcells.
#' There can be more detected flowcells than there are expected flowcells.
#' 
#' @param detected_flowcells A dataframe with the columns "flowcell_name" and "flowcell_lane".
#' @param expected_flowcells A dataframe with the columns "flowcell_name" and "flowcell_lane".
validate_detected_flowcells= function(detected_flowcells, expected_flowcells) {
  missing_flowcells= expected_flowcells %>% dplyr::anti_join(detected_flowcells, by= c('flowcell_name', 'flowcell_lane'))

  if(nrow(missing_flowcells) != 0) {
    print('WARNING: The following flowcells/lanes specified in the sample meta were not detected in the fastq reads.')
    print(missing_flowcells)
    print('Check that the sample meta is correct or that all fastq files are in the correct directory.')
  }
}

#' process_in_chunks
#' 
#' This function runs some action over chunks of a large file. At the end, all chunks are
#' appended together.
#' 
#' @param large_file_path description
#' @param chunk_size description
#' @param action A function passed to act on each chunk
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
  
  output_table= data.table::rbindlist(chunk_collector)
  return(output_table)
}

#' collate_fastq_reads
#' 
#' This function takes in the fastq reads (uncollapsed_raw_counts) and
#' filters for reads coming from flowcells specified in the sample meta.
#' The function then sums up the reads across specified sequencing index columns and
#' maps the sequencing index columns to the ID columns.
#' 
#' @param uncollapsed_raw_counts Dataframe of reads from all the fastq files with the following columns -
#'                    "flowcell_name", "flowcell_lane", "index_1", "index_2", "forward_read_cl_barcode", and "n".
#'                    The flowcell columns are optional. If they do not exists, flowcell filters will be skipped.
#' @param sample_meta Sample metadata generate for the project which may contain the following columns - 
#'                    "flowcell_names", "flowcell_lanes", "index_1", "index_2". The sample meta MUST contain
#'                    "flowcell_names" and "flowcell_lanes" for filtering.
#' @param sequencing_index_cols Sequencing columns from the sample meta that the counts should be collapsed on.
#'                              These columns should be a subset of the four sequencing related columns in the
#'                              sample meta - "flowcell_names", "flowcell_lanes", "index_1", and "index_2". They 
#'                              should also uniquely identify every PCR well. This parameter defaults onto 
#'                              the following columns: "index_1", "index_2".
#' @param id_cols ID columns from the sample meta that uniquely identify every PCR well. These columns should not 
#'                include any sequencing related columns. This parameter defaults onto "pcr_plate", "pcr_well". This 
#'                parameter can also be a list of the sample conditions columns as long as they uniquely identify every
#'                PCR well. For example "cell_set", "treatment", "dose", "day", "bio_rep", "tech_rep" can also be used.
#' @param reverse_index2 Index 2 should be reversed if the sequencer uses a reverse complement workflow. 
#'                       Defaults to FALSE.
#' @param barcode_col String name of the column in uncollapsed_raw_counts that contains the sequences.  
#' @returns Returns a dataframe with columns specified by the id_cols along with barcode_col, and "n".
#' @import tidyverse
#' @import data.table
collate_fastq_reads= function(uncollapsed_raw_counts, sample_meta, 
                              sequencing_index_cols= c('index_1', 'index_2'),
                              id_cols= c('pcr_plate', 'pcr_well'),
                              reverse_index2= FALSE,
                              barcode_col= 'forward_read_cl_barcode') {
  require(tidyverse)
  require(data.table)
  
  # Reverse index 2 if specified ----
  if(reverse_index2) {
    if('index_2' %in% colnames(sample_meta)) {
      print("Reverse-complementing index 2 barcode ...")
      sample_meta$index_2= chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$index_2))
    } else {
      stop('Reverse index 2 is set to TRUE, but index_2 does not exists.')
    }
  }
  
  # Create sequence map ----
  sequencing_map= sample_meta %>% dplyr::distinct(pick(all_of(c(sequencing_index_cols, id_cols))))
  
  # Validation: Check that flowcell_names and flowcell_lanes exist in the sample meta ----
  if(!validate_columns_exist(c('flowcell_names', 'flowcell_lanes'), sample_meta)) {
    stop('The above column(s) are NOT present in the sample meta.')
  }
  
  # Validation: Check that sequencing_index_cols exist in the sample meta ----
  if(!validate_columns_exist(sequencing_index_cols, sample_meta)) {
    print('The following sequencing_index_cols are not present in the sample meta.')
    stop('The above sequencing_index_cols are NOT present in the sample meta.')
  }
  
  # Validation: Check that id_cols exist in the sample meta ----
  if(!validate_columns_exist(id_cols, sample_meta)) {
    stop('One or more id_cols is NOT present in the sample meta.')
  }
  
  # Validation: Check that sequencing_index_cols in the sample meta are filled out ----
  # Check for rows in sequencing_index_cols that equate to empty - NA, "NA", "", " "
  # Error out of the sequencing_index_cols are not filled out in the sample meta.
  if(!validate_columns_entries(sequencing_index_cols, sample_meta)) {
    stop('One or more sequencing_index_cols in the sample meta is not filled out.')
  }
  
  # Validation: Check that mapping is one to one ----
  check_mapping= sequencing_map %>% dplyr::group_by(pick(all_of(sequencing_index_cols))) %>%
    dplyr::filter(dplyr::n() > 1) %>% dplyr::ungroup()
  if(nrow(check_mapping) > 0) {
    print('The following sequencing locations map to multiple conditions.')
    print(check_mapping)
    stop('The sequencing index columns do not map 1 to 1 to the ID columns.')
  }
  
  # If "flowcell_name" and "flowcell_lane" are present, filter for valid flowcells ----
  # Note: Can this switch be tied to the sequencer type?
  if(base::all(c('flowcell_name', 'flowcell_lane') %in% colnames(uncollapsed_raw_counts))) {
    print('Detecting flowcells. Filtering for valid flowcells ...')
    
    # Determine which flowcell names + lanes are expected ----
    # "flowcell_names" and "flowcell_lanes" are strings that can contain more than one item.
    # Columns can be parsed by splitting on the chars , ; :
    # If there are multiple lane names and lane numbers, this uses the Cartesian product!
    # Note: fread and read.csv keeps commas, read_csv DROPS commas
    expected_flowcells= sample_meta %>% dplyr::distinct(flowcell_names, flowcell_lanes) %>%
      dplyr::mutate(flowcell_name= base::strsplit(flowcell_names, split='[,;:]', fixed=F),
                    flowcell_lane= base::strsplit(flowcell_lanes, split='[,;:]', fixed=F)) %>% 
      tidyr::unnest(cols= flowcell_name) %>% tidyr::unnest(cols= flowcell_lane) %>%
      dplyr::mutate(flowcell_lane= as.numeric(flowcell_lane))
    
    # Print out expected flowcells from the sample meta.
    print(paste0('Identified ', nrow(expected_flowcells), ' unique flowcell + lane combos in the sample meta ...'))
    print(expected_flowcells)
    
    # Print warning if there are multiple flowcell names with multiple flowcell lanes.
    multi_name_and_lanes= expected_flowcells %>% dplyr::filter(grepl(',:;', flowcell_names) & grepl(',:;', flowcell_names))
    if(nrow(multi_name_and_lanes) > 0) {
      print('WARNING: Detected sample(s) sequenced over multiple flowcells and flowcell lanes.')
      print('The function assumes that the same lanes were used for both flowcells.')
    }
    
    # Validation: Check that all expected flowcell name + lanes are detected ----
    # Check that all expected flowcell name + lanes are present in uncollapsed raw counts.
    detected_flowcells= uncollapsed_raw_counts %>% dplyr::distinct(flowcell_name, flowcell_lane)
    print(paste0('Identified ', nrow(detected_flowcells), ' unique flowcell + lane combos in the uncollapsed raw counts ...'))
    print(detected_flowcells)
    validate_detected_flowcells(detected_flowcells, expected_flowcells)
    
    # Validation: Check that sequencing_index_cols uniquely identify rows of sample meta ----
    if(!validate_unique_samples(sequencing_index_cols, sample_meta)) {
      print('There may be multiple entries in the sample meta that have the same combination of sequencing index columns.')
      stop('The specified sequencing index columns do NOT uniquely identify every PCR well.')
    }
    
    # Validation: Check that id_cols uniquely identify rows of sample meta ----
    if(!validate_unique_samples(id_cols, sample_meta)) {
      print('There may be multiple entries in the sample meta that have the same combination of ID columns.')
      stop('The specified ID columns do NOT uniquely identify every PCR well.')
    }
    
    # Filter for expected flowcells ----
    uncollapsed_raw_counts= data.table::merge.data.table(
      uncollapsed_raw_counts, data.table::setDT(expected_flowcells), 
      by= c('flowcell_name', 'flowcell_lane'), allow.cartesian= FALSE)
    
  } else {
    print('Flowcell_name and/or flowcell_lane were not detected in raw_counts_uncollapsed.')
    print('Proceeding without filtering flowcells ...')
  }
  
  # Create summed_reads file ----
  # Filter for the expected flowcells and summed up the reads over the ID cols.
  print('Summing up reads ...')
  summed_reads= data.table::merge.data.table(uncollapsed_raw_counts, sequencing_map, by= sequencing_index_cols)
  summed_reads= summed_reads[, .(n= sum(n)), by= c(id_cols, barcode_col)]
  
  # Escape for when a chunk contains invalid sequencing locations
  if(nrow(summed_reads) == 0) {
    print('WARNING: summed_reads is empty!')
    return(summed_reads)
  }
  
  # Calculate index purity in a chunk----
  index_purity= sum(summed_reads$n) / sum(uncollapsed_raw_counts$n)
  print(paste0('Index purity in chunk: ', round(index_purity, 4)))
  if(index_purity > 1) {
    stop('ERROR: Chunk index purity is greater than 1!')
  } else if(index_purity < 0.5) {
    print('Warning: Low index purity!')
  } else {}
  
  print('Collate_fastq_reads has completed!')
  return(summed_reads)
}

#' extract_known_barcodes
#' 
#' This function runs some action over chunks of a large file. At the end, all chunks are
#' appended together.
#' 
#' @param raw_counts description
#' @param known_barcodes A vector known barcodes.
#' @param barcode_col String name of the column in uncollapsed_raw_counts that contains the sequences. 
extract_known_barcodes= function(summed_reads, known_barcodes, barcode_col= 'forward_read_cl_barcode') {
  # Create boolean column of known or unknown
  summed_reads[, known := get(barcode_col) %chin% known_barcodes]
  
  # Filter using that boolean column
  unknown_reads= summed_reads[known == FALSE,][order(-n)][, known := NULL]
  summed_reads= summed_reads[known == TRUE,][, known := NULL]
  
  return(list(unknown_reads= unknown_reads, known_reads= summed_reads))
}
