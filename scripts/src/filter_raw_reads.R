suppressPackageStartupMessages(library(sets))

#' filter raw reads
#' 
#' takes the raw readcount table and filters for expected indices and cell lines
#' using the given metadata. QC metrics are returned as a data.frame
#'
#' @param raw_counts - an unfiltered counts table
#' @param sample_meta - the sample metadata for the particular experiment. Must follow the given set of 
#'                      guidelines for metadata. Required columns include:
#'                      - IndexBarcode1
#'                      - IndexBarcode2
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
#' @param id_cols - vector of column names used to generate unique profile_id for each sample. 
#'                  cell_set,treatment,dose,dose_unit,day,bio_rep,tech_rep by default
#' @param reverse_index2 Reverses index2 for certain sequencers
#' @param count_threshold Threshold to call low counts. 
#' @return - list with the following elements
#' #' \itemize{
#'   \item unmapped_reads: table of reads with valid index pairs but did not map to any known barcode.
#'         The table contains the following columns - index_1, index_2, forward_read_cl_barcode, and n
#'   \item annotated_counts: table of reads and the associated well and well conditions
#'   \item filtered_counts: table of all expected reads for the project, this is a subset of annotated counts.
#'   \item qc_table: QC table of index_purity and cell_line_purity 
#' }
#' @export 
filter_raw_reads = function(
  raw_counts, sample_meta, cell_line_meta, 
  cell_set_meta, CB_meta,
  id_cols=c('cell_set','treatment','dose','dose_unit','day','bio_rep','tech_rep'),
  reverse_index2= FALSE, count_threshold= 40) {
  
  require(tidyverse)
  
  # Processing metadata ---- 
  print("Converting CB_meta from log10 to log2 ...")
  CB_meta= CB_meta %>% dplyr::mutate(log2_dose= log_dose/log10(2)) %>% dplyr::select(-log_dose)
  
  if (reverse_index2) {
    print("Reverse-complementing index 2 barcode ...")
    sample_meta$IndexBarcode2 <- chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$IndexBarcode2))
  }
  
  # Filtering by index barcodes ----
  print("Filtering raw counts ...")
  index_filtered= raw_counts %>% 
    dplyr::filter(index_1 %in% sample_meta$IndexBarcode1, index_2 %in% sample_meta$IndexBarcode2) %>%
    dplyr::mutate(mapped= ifelse(forward_read_cl_barcode %in% c(cell_line_meta$Sequence, CB_meta$Sequence), T, F))
  
  # index purity for QC table
  index_purity= sum(index_filtered$n)/ sum(raw_counts$n)
  
  # Mapping barcodes with valid indices ----
  print('Splitting off unmapped reads ...')
  # unmapped reads are defined as having valid indices but do not map to barcodes in PRISM
  unmapped_reads= index_filtered %>% dplyr::filter(mapped==F) %>% dplyr::select(-mapped) %>% 
    dplyr::arrange(dplyr::desc(n))
  
  # make template of expected reads
  print('Creating template of expected reads ...')
  template= sample_meta %>% tidyr::unite('profile_id', all_of(id_cols), sep=':', remove=F, na.rm=F) %>%
    dplyr::left_join(cell_set_meta, by= 'cell_set') %>%
    dplyr::mutate(members= ifelse(is.na(members), str_split(cell_set, ';'), str_split(members, ';'))) %>% 
    tidyr::unnest(cols=c(members)) %>%
    dplyr::left_join(cell_line_meta, by= dplyr::join_by('members'=='LUA')) # NEW
  
  # check for control barcodes and add them to the template
  if(any(unique(sample_meta$control_barcodes) %in% c('Y', 'T', T))) {
    cb_template= sample_meta %>% dplyr::filter(control_barcodes %in% c('Y', 'T', T)) %>%
      dplyr::mutate(joiner= 'temp') %>%
      dplyr::inner_join(CB_meta %>% dplyr::mutate(joiner= 'temp'), by='joiner', relationship= 'many-to-many') %>% 
      dplyr::select(-joiner)
    template= dplyr::bind_rows(template, cb_template)
  }
  
  # Annotating reads ----
  print("Annotating reads ...")
  annotated_counts= index_filtered %>% dplyr::filter(mapped) %>% # New: drop unmapped reads
    dplyr::left_join(cell_line_meta, by= join_by('forward_read_cl_barcode'=='Sequence')) %>%
    dplyr::left_join(CB_meta, by= join_by('forward_read_cl_barcode'=='Sequence')) %>%
    dplyr::left_join(sample_meta, by= join_by('index_1'=='IndexBarcode1', 'index_2'=='IndexBarcode2')) %>%
    merge(template %>% dplyr::mutate(expected_read= T),
          by.x= c('index_1', 'index_2', 'forward_read_cl_barcode', intersect(colnames(template), colnames(.))), 
          by.y= c('IndexBarcode1', 'IndexBarcode2', 'Sequence', intersect(colnames(template), colnames(.))), 
          all.x=T, all.y=T) %>% 
    # drop unneeded columns and fill in any new NAs from the merge
    dplyr::select(!any_of(c('prism_cell_set', 'members', 'mapped'))) %>%
    dplyr::mutate(n= replace_na(n, 0), expected_read= replace_na(expected_read, F))
  
  # Generating filtered reads ----
  print("Filtering reads ...")
  filtered_counts= annotated_counts %>% dplyr::filter(expected_read) %>%
    dplyr::select(!any_of(c('index_1', 'index_2', 'forward_read_cl_barcode', 'LUA', 'expected_read'))) %>%
    dplyr::mutate(flag= ifelse(n==0, 'Missing', NA),
                  flag= ifelse(n!=0 & n < count_threshold, 'low counts', flag))
  
  # cell line purity for QC table
  cell_line_purity = sum(filtered_counts$n)/ sum(index_filtered$n)
  
  # Generating QC table ----
  qc_table = data.frame(index_purity= index_purity, cell_line_purity= cell_line_purity)
  
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



