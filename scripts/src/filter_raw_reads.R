options(cli.unicode = FALSE)

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
#' Takes the raw readcount table and filters for expected indices and cell lines
#' using the given metadata.
#'
#' @param raw_counts Dataframe of reads. The columns of this dataframe should include the id_cols,
#'                   "forward_read_cl_barcode", and "n".
#' @param sample_meta Dataframe of the metadata for the sequencing run. This file should contain the id_cols,
#'                    "cell_set", "control_barcodes", etc.
#' @param cell_set_meta Master metadata of cell sets and their contents with the following required columns -
#'                      "cell_set" and "members".
#' @param cell_line_meta Master metadata of cell lines with the following required columns - "CCLE_name",
#'                       "DepMap_ID", "LUA", and "Sequence".
#' @param CB_meta Master metadata of control barcodes, their sequences, and their doses. The file should contain 
#'                the columns - "Sequence", "cb_name", and "cb_log10_dose".
#' @param id_cols Columns present in both raw_counts and sample_meta that uniquely identify each PCR well. 
#'                This defaults to "pcr_plate", "pcr_well".
#' @returns List with the following elements:
#' #' \itemize{
#'   \item annotated_counts: table of reads and the associated well and well conditions.
#'   \item filtered_counts: table of all expected reads for the project, this is a subset of annotated counts.
#' }
filter_raw_reads= function(prism_barcode_counts, 
                            sample_meta, cell_set_meta, cell_line_meta, CB_meta,
                            id_cols= c('pcr_plate', 'pcr_well'),
                            barcode_col= 'forward_read_cl_barcode') {
  require(magrittr)
  require(tidyverse)
  require(data.table)
  
  # CB meta may be in log10 and should be converted to log2 ----
  # Check if 'cb_log10_dose' or 'cb_log2_dose' is present in CB_meta.
  # If 'cb_log10_dose' is present, convert it to 'cb_log2_dose'.
  # If cb_log2_dose is present, print out an acknowledgement
  # If neither names are present, error out
  if('cb_log10_dose' %in% colnames(CB_meta)) {
    print("Converting CB_meta from log10 to log2.")
    CB_meta= CB_meta %>% dplyr::mutate(cb_log2_dose= cb_log10_dose/log10(2)) %>% dplyr::select(-cb_log10_dose)
  } else if('cb_log10_dose' %in% colnames(CB_meta)) {
    print('Detecting cb_log2_dose in CB_meta.')
  } else {
    stop('Missing either cb_log2_dose or cb_log10_dose in CB_meta.csv.')
  }
  
  # Validation: Check that id_cols exist in the sample meta ----
  if(!validate_columns_exist(id_cols, sample_meta)) {
    stop('One or more id_cols is NOT present in the sample meta.')
  }
  
  # Validation: Check that id_cols uniquely identify every rows of sample meta ----
  if(!validate_unique_samples(id_cols, sample_meta)) {
    print('There may be multiple entries in the sample meta that have the same combination of id_cols.')
    stop('The specified ID columns do NOT uniquely identify every PCR well.')
  }
  
  # Validation: Check that cell sets do not contain duplicate LUAs ----
  # This will produce a warning if a LUA appears in a cell set more than once!
  # This currently does NOT result in an error. Error avoided using a distinct when creating the template of expected reads.
  validate_cell_set_luas(sample_meta, cell_set_meta)
  
  # Creating a template of all expected reads in the run ----
  # Use all 4 metadata files to create a "template" dataframe where
  # every row is a cell line that is expected in a PCR well. 
  print('Creating template of expected reads.')
  
  # From the sample meta, identify all expected cell line sequences
  # Filter for just wells that have a specified cell set and join cell_set_and_pool_meta
  template= data.table::merge.data.table(sample_meta[!cell_set %in% c(NA, 'NA', '', ' '),],
                                         cell_set_and_pool_meta, by= 'cell_set', allow.cartesian= TRUE)
  # Left join barcode sequence using data.table inplace merge
  template[cell_line_meta, c(barcode_col) := get(barcode_col), on= 'depmap_id']
  
  # Check for control barcodes and add them to the template.
  if(any(unique(sample_meta$control_barcodes) %in% c('Y', 'T', T))) {
    # From the sample meta, identify all expected control barcode sequences
    # Filter for just wells that have cb_ladder(s) present in CB_meta and join the CB_meta
    cb_template= data.table::merge.data.table(sample_meta[cb_ladder %in% unique(CB_meta$cb_ladder),],
                                              CB_meta, by= 'cb_ladder', allow.cartesian= TRUE)
    template= data.table::rbindlist(list(template, cb_template), fill= TRUE)
  }
  
  # Add a column to indicate if a read was expected - will help in annotated counts
  template[, expected_read := TRUE]
  
  # Annotating reads ----
  # From prism_barcode_counts, left join metadata to annotate all reads.
  # Perform a full join with the template of expected reads so that there is a row entry for 
  # cell lines not detected in sequencing.
  print("Annotating reads.")
  # Annotate prism_barcode_counts using data.table left joins performed in place!
  # Left join sample meta using id_cols as the keys
  # Sample meta must be the first to join! Otherwise you risk cb_ladder and pool_id complicating join with the template.
  prism_barcode_counts[sample_meta, base::setdiff(colnames(sample_meta), id_cols) :=
                         mget(base::setdiff(colnames(sample_meta), id_cols)), on= id_cols]
  # Left join cell_line_meta using barcode_col as the key
  prism_barcode_counts[cell_line_meta, base::setdiff(colnames(cell_line_meta), barcode_col) :=
                         mget(base::setdiff(colnames(cell_line_meta), barcode_col)), on= barcode_col]
  # Left join CB_meta using barcode_col as the key
  prism_barcode_counts[CB_meta, base::setdiff(colnames(CB_meta), barcode_col) :=
                         mget(base::setdiff(colnames(CB_meta), barcode_col)), on= barcode_col]

  # Create annotated counts by performing a full join with template
  annotated_counts= data.table::merge.data.table(prism_barcode_counts, template,
                                                 by= intersect(colnames(prism_barcode_counts), colnames(template)),
                                                 all.x= TRUE, all.y= TRUE, allow.cartesian= FALSE)

  # Use data.table inplace compute fast ifelse.
  annotated_counts[, n := data.table::fifelse(is.na(n), 0, n)] # Imput zeros for undetected barcodes
  # Fill in the expected_read column with FALSE
  annotated_counts[, expected_read := data.table::fifelse(is.na(expected_read), F, expected_read)] 
  
  # Generating filtered reads ----
  # Get filtered counts from annotated counts and drop a few select columns.
  print("Filtering reads ...")
  filtered_counts= annotated_counts[expected_read, ] %>%
    dplyr::select(!any_of(c('flowcell_names', 'flowcell_lanes', 'index_1', 'index_2', 
                            'forward_read_cl_barcode', 'expected_read')))
  
  # Calculate cell line purity ----
  cell_line_purity= sum(filtered_counts$n)/ sum(prism_barcode_counts$n)
  print(paste0('Cell line purity: ', round(cell_line_purity, 4)))
  if(cell_line_purity > 1) {
    stop('ERROR: Cell line purity is greater than 1!')
  }
  if(cell_line_purity < 0.5) {
    print('Warning: Low cell line purity!')
  }
  
  # Return both annotated_counts and filtered_counts ----
  print('Filter_raw_reads has completed!')
  return(list(annotated_counts= annotated_counts, filtered_counts= filtered_counts))
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

