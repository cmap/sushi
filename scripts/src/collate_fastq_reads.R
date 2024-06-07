#' collate_fastq_reads
#' 
#' This function takes in the fastq reads (uncollapsed_raw_counts) and
#' filters for reads coming from flowcells specificed in the sample meta.
#' The function then sums up the reads the flowcells using the index pairs.
#' 
#' @param uncollapsed_raw_counts Data frame of reads from all the fastq files with the following columns - \cr
#'                    "flowcell_name", "flowcell_lane", "index_1", "index_2", and "forward_read_cl_barcode", "n"
#' @param sample_meta Sample metadata generate for the project  
#' @param seq_cols Sequencing columns from the sample meta taht the counts should be collapsed on. \cr
#'                 This defaults onto the following columns: "IndexBarcode1", "IndexBarcode2"
#' @returns Returns a dataframe with the following columns - "index_1", "index_2", and "n"
#' @import tidyverse
collate_fastq_reads= function(uncollapsed_raw_counts, sample_meta, seq_cols= c('IndexBarcode1', 'IndexBarcode2')) {
  require(tidyverse)
  
  # Change flowcell_lane to single lane in uncollapsed_raw_counts ----
  # Uncollapsed raw counts's flowcell_lane is a single value, but 
  # sample meta's flowcell_lane is a string of numbers to parse.
  # These should be different columns in the event that you need to collapse over different groups of lanes.
  if(!'single_lane' %in% colnames(uncollapsed_raw_counts)) {
    uncollapsed_raw_counts= dplyr::rename(uncollapsed_raw_counts, single_lane= flowcell_lane)
  }
  
  # Determine which flowcell names + lanes are expected ----
  # "flowcell_lane" is a string containing a group of lanes, to be parsed by splitting on the chars , ; :
  # "single_lane" contains the individual lane
  # Note: fread and read.csv keeps commas, read_csv DROPS commas
  meta_flowcells= sample_meta %>% dplyr::distinct(flowcell_name, flowcell_lane) %>%
    dplyr::mutate(single_lane= base::strsplit(flowcell_lane, split='[,;:]', fixed=F)) %>%
    tidyr::unnest(cols= single_lane) %>%
    dplyr::mutate(single_lane= as.numeric(single_lane))
  
  # Print out expected flowcells from the sample meta.
  print(paste0('Identified ', nrow(meta_flowcells), ' unique flowcell + lane combos in the sample meta ...'))
  print(meta_flowcells)
  
  # QC: Check that all expected flowcell name + lanes are detected ----
  # Check that all expected flowcell name + lanes are present in uncollapsed raw counts.
  missing_flowcells= meta_flowcells %>%
    dplyr::anti_join(uncollapsed_raw_counts %>% dplyr::distinct(flowcell_name, single_lane),
                     by= c('flowcell_name', 'single_lane'))
  if(nrow(missing_flowcells)!= 0) {
    print('The following flowcells/lanes specified in the sample meta were not detected in the fastq reads.')
    print(missing_flowcells)
    print('Check that the sample meta is correct or that all fastq files are in the correct directory.')
    stop('ERROR: One or more flowcell specified in the sample meta was not detected.')
  }
  
  # QC: Check that seq_cols are present in the sample meta ----
  if(any(!seq_cols %in% colnames(sample_meta))) {
    stop('ERROR: One or more seq_col is not present in the sample meta.')
  }
  
  # QC: Check that seq_cols in the sample meta are filled out ----
  # Check for rows in seq_cols that equate to empty - NA, "NA", "", " "
  # Error out of the seq_cols are not filled out in the sample meta.
  missing_rows_seq_cols= sample_meta %>% dplyr::filter(if_any(all_of(seq_cols), ~ . %in% c(NA, 'NA', '', ' ')))
  if(nrow(missing_rows_seq_cols) > 0) {
    print('The following rows in the sample meta are not filled out for all of the sequencing columns.')
    print(missing_rows_seq_cols) # show the empty rows
    stop('ERROR: One or more seq col in the sample meta is not filled out.')
  }
  
  # QC: Check that seq_cols uniquely identify rows of sample meta ----
  # Find the unique combinations of seq_cols and check that it matches the number of row in the sample meta.
  unique_seq_col_vals= sample_meta %>% dplyr::distinct(pick(all_of(seq_cols)))
  if(nrow(unique_seq_col_vals) != nrow(sample_meta)) {
    print('There may be multiple entries in the sample meta that have the same combination of sequencing columns.')
    stop('ERROR: The specified sequencing columns do NOT uniquely identify every PCR well.')
  }
  
  # Convert seq_cols to equivalent names in uncollapsed file ----
  # Some columns that represent the same things are called different names 
  # between the sample meta and uncollapsed raw counts. The following code converts one name to the other. 
  converted_seq_cols= c() 
  for(item in seq_cols) {
    if(item=='IndexBarcode1') {
      converted_seq_cols= c(converted_seq_cols, 'index_1')
    } else if(item=='IndexBarcode2') {
      converted_seq_cols= c(converted_seq_cols, 'index_2')
    } else {
      converted_seq_cols= c(converted_seq_cols, item)
    }
  }
  
  # Create raw counts by summing over appropriate seq_cols ----
  # Use an inner join to collect reads with valid flowcell name/lane combinations, 
  # then sum reads across sequencing columns
  raw_counts= uncollapsed_raw_counts %>% 
    dplyr::inner_join(meta_flowcells, by= c('flowcell_name', 'single_lane'), relationship= 'many-to-one') %>%
    dplyr::group_by(pick(all_of(c(converted_seq_cols, 'forward_read_cl_barcode')))) %>% 
    dplyr::summarize(n= sum(n)) %>% dplyr::ungroup()
  
  return(raw_counts)
}
