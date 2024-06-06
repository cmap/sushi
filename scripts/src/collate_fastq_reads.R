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
collate_fastq_reads= function(uncollapsed_raw_counts, sample_meta, seq_cols) {
  require(tidyverse)
  
  # Determine which flowcell names + lanes are expected ----
  # Parse flowcell_lane by splitting on the chars , ; :
  # Note: fread and read.csv keeps commas, read_csv DROPS commas
  meta_flowcells= sample_meta %>% dplyr::distinct(flowcell_name, flowcell_lane) %>%
    dplyr::mutate(flowcell_lane= base::strsplit(flowcell_lane, split='[,;:]', fixed=F)) %>%
    tidyr::unnest(cols= flowcell_lane) %>%
    dplyr::mutate(flowcell_lane= as.numeric(flowcell_lane))
  
  print(paste0('Identified ', nrow(meta_flowcells), ' unique flowcell + lane combos in the sample meta ...'))
  print(meta_flowcells)
  
  # QC: Check that all expected flowcell name + lanes are detected ----
  missing_flowcells= uncollapsed_raw_counts %>% dplyr::distinct(flowcell_name, flowcell_lane) %>% dplyr::mutate(present=T) %>%
    dplyr::right_join(meta_flowcells, by=c('flowcell_name', 'flowcell_lane')) %>%
    dplyr::filter(is.na(present))
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
  # check for rows in seq_cols that equate to empty - NA, "NA", "", " "
  missing_rows_seq_cols= sample_meta %>% dplyr::filter(if_any(all_of(seq_cols), ~ . %in% c(NA, 'NA', '', ' ')))
  if(nrow(missing_rows_seq_cols) > 0) {
    print('The following rows in the sample meta are not filled out for all of the sequencing columns.')
    print(missing_rows_seq_cols)
    stop('ERROR: One or more seq col in the sample meta is not filled out.')
  }
  
  # Convert seq_cols to equivalent names in uncollapsed file ----
  # column names in uncollapsed_raw_counts are not the same as those in the sample meta
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
  # use an inner join to collect reads with valid flowcell name/lane combinations, then sum using index pairs
  raw_counts= uncollapsed_raw_counts %>% 
    dplyr::inner_join(meta_flowcells, by= c('flowcell_name', 'flowcell_lane'), relationship= 'many-to-one') %>%
    dplyr::group_by(pick(all_of(c(converted_seq_cols, 'forward_read_cl_barcode')))) %>% 
    dplyr::summarize(n= sum(n)) %>% dplyr::ungroup()
  
  return(raw_counts)
}
