#' collate_fastq_reads
#' 
#' This function takes in the fastq reads (uncollapsed_raw_counts) and
#' filters for reads coming from flowcells specificed in the sample meta.
#' The function then sums up the reads the flowcells using the index pairs.
#' 
#' @param sample_meta Sample metadata generate for the project
#' @param fastq_reads Data frame of reads from all the fastq files with the following columns - \cr
#'                    "flowcell_name", "flowcell_lane", "index_1", "index_2", and "n"
#' @returns Returns a dataframe with the following columns - "index_1", "index_2", and "n"
collate_fastq_reads= function(sample_meta, fastq_reads) {
  # determine which flowcell names + lanes are expected
  meta_flowcells= sample_meta %>% dplyr::distinct(flowcell_name, flowcell_lane) %>%
    tidyr::separate_longer_delim(flowcell_lane, delim=',') %>% # fread and read.csv keeps commas, read_csv DROPS commas
    dplyr::mutate(flowcell_lane= as.numeric(flowcell_lane))
  
  print(paste0('Identified ', nrow(meta_flowcells), ' unique flowcell + lane combos in the sample meta ...'))
  print(meta_flowcells)
  
  # check that expected flowcells are detected
  missing_flowcells= fastq_reads %>% dplyr::distinct(flowcell_name, flowcell_lane) %>% dplyr::mutate(present=T) %>%
    dplyr::right_join(meta_flowcells, by=c('flowcell_name', 'flowcell_lane')) %>%
    dplyr::filter(is.na(present))
  if(nrow(missing_flowcells)!= 0) {
    print('The following flowcells/lanes specified in the sample meta were not detected in the fastq reads.')
    print(missing_flowcells)
    print('Check that the sample meta is correct or that all fastq files are in the correct directory.')
    stop()
  }
  
  # use an inner join to collect reads with valid flowcell name/lane combinations, then sum using index pairs
  raw_counts= fastq_reads %>% dplyr::inner_join(meta_flowcells, by=c('flowcell_name', 'flowcell_lane')) %>%
    dplyr::group_by(index_1, index_2, forward_read_cl_barcode) %>% 
    dplyr::summarize(n= sum(n)) %>% dplyr::ungroup()
  
  return(raw_counts)
}
