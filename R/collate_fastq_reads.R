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
  # expects flowcell_lanes to be a string of numbers
  unique_flowcells= sample_meta %>% dplyr::distinct(flowcell_name, flowcell_lane) %>%
    tidyr::separate_longer_position(flowcell_lane, width=1) # read_csv and read.csv does NOT put commas between numbers!
  print(paste0('Identified ', nrow(unique_flowcells), ' unique flowcell + lane combos in the sample meta ...'))
  print(unique_flowcells)
  
  # use an inner join to collect reads with valid flowcell name/lane combinations, then sum using index pairs
  # merge might take a long time ...
  raw_counts= fastq_reads %>% base::merge(unique_flowcells, by=c('flowcell_name', 'flowcell_lane')) %>%
    dplyr::group_by(index_1, index_2, forward_read_cl_barcode) %>% 
    dplyr::summarize(n= sum(n)) %>% dplyr::ungroup()
  
  return(raw_counts)
}