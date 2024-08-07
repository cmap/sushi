options(cli.unicode = FALSE)

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

#' validate_unique_samples
#' 
#' This function checks that a list of columns uniquely identifies all entries of a dataframe.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_unique_samples= function(selected_columns, df) {
  message= paste0('The following columns do not uniquely identify every row of the dataframe: ',
                  paste(selected_columns, collapse=', '))
  print(message)
  unique_column_values= df %>% dplyr::distinct(pick(all_of(selected_columns)))
  if(nrow(unique_column_values) != nrow(df)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' validate_cell_set_luas
#' 
#' This function checks that every cell set in the sample meta does not contain duplicate members.
#' If a cell set has duplicate LUAs, a warning is printed.
#' 
#' @param sample_meta The sample_meta df with the column "cell_set".
#' @param cell_set_meta The cell_set_meta df with the columns "cell_set" and "members".
validate_cell_set_luas= function(sample_meta, cell_set_meta) {
  duplicate_luas= cell_set_meta %>% dplyr::filter(cell_set %in% unique(sample_meta$cell_set)) %>%
    dplyr::mutate(members= str_split(members, ';')) %>%
    tidyr::unnest(cols= 'members') %>%
    dplyr::count(cell_set, members, name= 'count') %>% dplyr::filter(count > 1)
  
  if(nrow(duplicate_luas) > 0) {
    print('WARNING: The following LUAs appear more than once in a cell set!!!')
    print(duplicate_luas)
    print('The module will continue, but check the cell_set meta!!!')
  }
}

#' filter raw reads
#' 
#' takes the raw readcount table and filters for expected indices and cell lines
#' using the given metadata. QC metrics are returned as a data.frame
#'
#' @param raw_counts - an unfiltered counts table
#' @param sample_meta - the sample metadata for the particular experiment. Must follow the given set of 
#'                      guidelines for metadata. Required columns include:
#'                      - index_1
#'                      - index_2
#'                      - cell_set
#' @param cell_line_meta - master metadata of cell lines with the following required columns:
#'                       - CCLE_name
#'                       - DepMap_ID
#'                       - LUA
#'                       - Sequence
#' @param cell_set_meta - master metdata of cell sets and their contents with the following required columns:
#'                      - cell_set
#'                      - members
#' @param CB_meta - master metdata of control barcodes, their sequences, and their doses. 
#'                The file should contain the columns:
#'                - Sequence
#'                - Name
#'                - log_dose
#' @param sequencing_index_cols Vector of column names from the sample meta that is used to uniquely identify where a
#'                              sequencing read is coming from. This defaults to "index_1" and "index_2", but 
#'                              it can be expanded to include "flowcell_names" and "flowcell_lanes".
#' @param reverse_index2 Reverses index2 for certain sequencers
#' @param count_threshold Threshold to call low counts. 
#' @param control_type - how the negative controls are designated in the trt_type column in the sample metadata
#' @return - list with the following elements
#' #' \itemize{
#'   \item unmapped_reads: table of reads with valid index pairs but did not map to any known barcode.
#'         The table contains the following columns - index_1, index_2, forward_read_cl_barcode, and n
#'   \item annotated_counts: table of reads and the associated well and well conditions
#'   \item filtered_counts: table of all expected reads for the project, this is a subset of annotated counts.
#'   \item qc_table: QC table of index_purity and cell_line_purity 
#' }
#' @export 
filter_raw_reads = function(raw_counts, 
                            sample_meta, cell_line_meta, cell_set_meta, CB_meta,
                            reverse_index2,
                            sequencing_index_cols= c('index_1', 'index_2'),
                            count_threshold= 40) {
  require(tidyverse)
  require(magrittr)
  
  # Processing metadata and inputs ---- 
  # CB meta is in log10 and should be converted to log2.
  print("Converting CB_meta from log10 to log2 ...")
  CB_meta= CB_meta %>% dplyr::mutate(log2_dose= log_dose/log10(2)) %>% dplyr::select(-log_dose)
  
  if (reverse_index2) {
    if ('index_2' %in% colnames(sample_meta)) {
      print("Reverse-complementing index 2 barcode ...")
      sample_meta$index_2 <- chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$index_2))
    } else {
      stop('ERROR: Reverse index 2 is set to TRUE, but index_2 does not exists.')
    }
  }
  
  # Validation: Check that sequencing_index_cols exist in the sample meta ----
  if(!validate_columns_exist(sequencing_index_cols, sample_meta)) {
    stop('One or more sequencing_index_cols is NOT present in the sample meta.')
  }
  
  # Validation: Check that sequencing_index_cols uniquely identify every rows of sample meta ----
  if(!validate_unique_samples(sequencing_index_cols, sample_meta)) {
    print('There may be multiple entries in the sample meta that have the same combination of sequencing index columns.')
    stop('The specified sequencing index columns do NOT uniquely identify every PCR well.')
  }
  
  # Validation: Check that cell sets do not contain duplicate LUAs ----
  # This will produce a warning if a LUA appears in a cell set more than once!
  # This currently does NOT result in an error. Error avoided using a distinct later in line 162
  validate_cell_set_luas(sample_meta, cell_set_meta)
  
  # Filtering by sequencing columns ----
  # Filter raw counts using the sequencing columns. 
  # Also create "mapped" column to identify reads that mapped to all known PRISM sequences.
  print("Filtering by sequencing columns ...")
  unique_sequencing_index_vals= sample_meta %>% dplyr::distinct(pick(all_of(sequencing_index_cols)))
  index_filtered= raw_counts %>% dplyr::semi_join(unique_sequencing_index_vals, by= sequencing_index_cols) %>%
    dplyr::mutate(mapped= forward_read_cl_barcode %in% c(cell_line_meta$Sequence, CB_meta$Sequence))
  
  # Calculate index purity for QC table.
  index_purity= sum(index_filtered$n)/ sum(raw_counts$n)
  
  # Split off unmapped reads ----
  # Unmapped reads are defined as having valid indices but do not map to barcodes in PRISM.
  # Also sorted reads in descending order by read count.
  print('Splitting off unmapped reads ...')
  unmapped_reads= index_filtered %>% dplyr::filter(mapped==F) %>% dplyr::select(-mapped) %>% 
    dplyr::arrange(dplyr::desc(n))
  
  # Creating a template of all expected reads in the run ----
  # Use all 4 meta data files to create a "template" dataframe where
  # every row is a cell line that is expected in a PCR well. 
  print('Creating template of expected reads ...')
  # Join cell_set_meta and cell_line_meta. The cell_set can be a name "P939" or a list of LUAs.
  template= sample_meta %>% dplyr::left_join(cell_set_meta, by= 'cell_set') %>%
    dplyr::mutate(members= ifelse(is.na(members), str_split(cell_set, ';'), str_split(members, ';'))) %>% 
    tidyr::unnest(cols= members) %>%
    dplyr::left_join(cell_line_meta, by= dplyr::join_by('members'=='LUA'), relationship= 'many-to-one') %>%
    dplyr::distinct() # To remove duplicate cell lines - cell sets sometimes have the same cell line repeated.
  # May need to replace distinct with an error later!
  
  # Check for control barcodes and add them to the template.
  if(any(unique(sample_meta$control_barcodes) %in% c('Y', 'T', T))) {
    cb_template= sample_meta %>% dplyr::filter(control_barcodes %in% c('Y', 'T', T)) %>%
      dplyr::mutate(joiner= 'temp') %>%
      dplyr::inner_join(CB_meta %>% dplyr::mutate(joiner= 'temp'), by='joiner', relationship= 'many-to-many') %>% 
      dplyr::select(-joiner)
    template= dplyr::bind_rows(template, cb_template)
  }
  
  # Annotating reads ----
  # From the set of reads that have the valid sequencing_index_cols combinations and map to the PRISM seq library,
  # join in metadata to give each read a name and PCR location.
  # Reads that to not match to the template are contaminants and,
  # reads that are only present in the template are missing/not detected by PCR.
  print("Annotating reads ...")
  annotated_counts= index_filtered %>% dplyr::filter(mapped) %>%
    dplyr::left_join(cell_line_meta, by= join_by('forward_read_cl_barcode'=='Sequence'),
                     relationship= 'many-to-one') %>%
    dplyr::left_join(CB_meta, by= join_by('forward_read_cl_barcode'=='Sequence'),
                     relationship= 'many-to-one') %>%
    dplyr::left_join(sample_meta, by= sequencing_index_cols, relationship= 'many-to-one') %>%
    dplyr::full_join(template %>% dplyr::mutate(expected_read= T),
                     by= c('forward_read_cl_barcode'='Sequence', intersect(colnames(template), colnames(.))),
                     relationship= 'one-to-one') %>%
    # drop unneeded columns and fill in any new NAs from the merge
    dplyr::select(!any_of(c('prism_cell_set', 'members', 'mapped'))) %>%
    dplyr::mutate(n= replace_na(n, 0), expected_read= replace_na(expected_read, F))
  
  # Generating filtered reads ----
  # Get filtered counts from annotated counts. Also flag reads that are either missing,
  # or below a count threshold.
  print("Filtering reads ...")
  filtered_counts= annotated_counts %>% dplyr::filter(expected_read) %>%
    dplyr::select(!any_of(c('flowcell_names', 'flowcell_lanes', 'index_1', 'index_2', 
                            'forward_read_cl_barcode', 'LUA', 'expected_read'))) %>%
    dplyr::mutate(flag= ifelse(n==0, 'Missing', NA),
                  flag= ifelse(n!=0 & n < count_threshold, 'low counts', flag))
  
  # Calculate cell line purity for the QC table.
  cell_line_purity= sum(filtered_counts$n)/ sum(index_filtered$n)
  
  # Generating QC table ----
  print('Generating QC table ...')
  qc_table= data.frame(index_purity= index_purity, cell_line_purity= cell_line_purity)
  
  return(list(unmapped_reads= unmapped_reads, annotated_counts= annotated_counts, 
              filtered_counts= filtered_counts, qc_table= qc_table))
}

# checks is a string can be numeric
is_numeric_string <- function(string){
  # tries converting the string to a number
  numeric_value <- suppressWarnings(as.numeric(string))
  
  # if the conversion was successful and if the input is not NA
  if (!sum(is.na(numeric_value))) {
    return(TRUE)  # String can be converted to a number
  } else {
    return(FALSE) # String cannot be converted to a number
  }
}

# checks if a string can be numeric and then converts it to numeric
convert_string_num_to_numeric <- function(df){
  # if a string can be a number, make it a number because it will be a number in filtered counts
  to_num <- sapply(df, is_numeric_string) # checks if can be numeric
  if(sum(to_num) != 0){ # make into number
    can_number <- which(to_num == T)
    for(ind in can_number){
      df[,ind] <- as.numeric(df[,ind])
    }
  }
  return(df)
}

remove_data = function(filtered_counts, data_to_remove) {
  filt_rm <- filtered_counts
  num_col_rm <- length(colnames(data_to_remove))
  
  # remove data_to_remove has any NAs, NULL or empty spaces, throw error
  if(sum(is.na(data_to_remove)) != 0){
    print("ERROR: NAs in removal data. Please fix this and try again.")
  }else if(any(sapply(data_to_remove, function(x) any(x == ""))) | any(sapply(data_to_remove, function(x) any(x == " ")))){
    print("ERROR: Empty spaces in removal data. Please fix this and try again.")
  }else if(any(sapply(data_to_remove, is.null))){
    print("ERROR: NULL spaces in removal data. Please fix this and try again.")
  }
  
  # make sure column names are in filtered_counts
  if(sum(colnames(data_to_remove) %in% colnames(filtered_counts)) != num_col_rm){
    print("ERROR: data_to_remove.csv columns names are not in filtered_counts.csv")
  }
  
  print("Removing data...")
  
  # split data_to_remove into two df: wild_df (with wild string), full_df (no wild string in any row)
  wild_string <- "EVERY"
  wild_df <- data_to_remove %>%
    filter(rowSums(sapply(., function(x) grepl(wild_string, x))) > 0)
  full_df <- anti_join(data_to_remove, wild_df)
  
  # anti_join the full data
  if(nrow(full_df) != 0){
    # if a string can be a number, make it a number because it will be a number in filtered counts
    full_df <- convert_string_num_to_numeric(full_df)
    filt_rm <- anti_join(filt_rm, full_df)
  }
  
  # anti_join the data with the wild string
  if(nrow(wild_df) != 0){
    # row-wise anti_join
    for(r in 1:dim(wild_df)[1]){
      crit_to_rm <- wild_df[r,]
      #print(crit_to_rm)
      
      # if there is an "EVERY", remove that column
      crit_to_rm <- crit_to_rm %>% 
        select(which(colSums(crit_to_rm == wild_string) == 0))
      
      # if a string can be a number, make it a number because it will be a number in filtered counts
      crit_to_rm <- convert_string_num_to_numeric(crit_to_rm)
      filt_rm <- anti_join(filt_rm, crit_to_rm) 
    } 
  }
  
  return(filt_rm)
}

