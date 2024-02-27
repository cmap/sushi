#' write_df_from_fastq
#' 
#' takes PRISM miseq or hiseq fastqs and returns a dataframe
#' with the number of barcodes from each cell line in each well
#' 
#' @param forward_read_fastq_files - vector of fastq file paths
#' @param index_1_files - vector of fastq file paths
#' @param index_2_files - vector of fastq file paths
#' @param write_interval - integer for how often a temp count file is written, NA by default. 
#' @return - cumulative_count_df A data.frame of readcounts by index_1, index_2 and forward_read_cl_barcode
#' @export 
write_df_from_fastq <- function(
  forward_read_fastq_files, 
  index_1_files, 
  index_2_files,
  write_interval = NA,
  CL_BC_LENGTH = 24,
  PLATE_BC_LENGTH = 8,
  WELL_BC_LENGTH = 8,
  save_loc = getwd()) 
{
  write_interval = as.numeric(write_interval)
  
  n_total_reads <- 0
  #cumulative_count_df <- data.frame()
  cumulative_count_df <- list() # new

  
  lp = 1 # new
  for (i in 1:length(forward_read_fastq_files)) {
    forward_stream <- ShortRead::FastqStreamer(forward_read_fastq_files[i])
    index_1_stream <- ShortRead::FastqStreamer(index_1_files[i])
    index_2_stream <- ShortRead::FastqStreamer(index_2_files[i])
    
    print(paste0("forward read file is: ", forward_read_fastq_files[i]))
    print(paste0("index 1 file is: ", index_1_files[i]))
    print(paste0("index 2 file is: ", index_2_files[i]))
    
    j <- 0
    repeat{
      forward_reads_chunk <- ShortRead::yield(forward_stream)
      index_1_chunk <- ShortRead::yield(index_1_stream)
      index_2_chunk <- ShortRead::yield(index_2_stream)
      
      if (length(forward_reads_chunk) == 0) {break}
      
      print(paste0('Processing read ', j * 1e6, ' to ', 
                   j*1e6 + length(forward_reads_chunk), 
                   ' from file ', i)) #(ShortRead)The default size for both streams and samples is 1M records;
      j <- j + 1
      
      forward_reads_string_set <- forward_reads_chunk  %>%
        ShortRead::sread() %>%
        Biostrings::DNAStringSet()
      
      index_1_string_set <- index_1_chunk %>%
        ShortRead::sread() %>%
        Biostrings::DNAStringSet()
      index_2_string_set <- index_2_chunk %>%
        ShortRead::sread() %>%
        Biostrings::DNAStringSet()
      
      forward_reads_cl_barcode <- XVector::subseq(forward_reads_string_set, 1, CL_BC_LENGTH)
 
      index_1 <- XVector::subseq(index_1_string_set, start = 1, end = PLATE_BC_LENGTH)
      index_2 <- XVector::subseq(index_2_string_set, start = 1, end = WELL_BC_LENGTH)
      
      forward_reads_cl_barcode = as.character(forward_reads_cl_barcode)
      index_1 = as.character(index_1)
      index_2 = as.character(index_2)
      chunk_df <- data.frame(forward_read_cl_barcode = forward_reads_cl_barcode, index_1 = index_1, index_2 = index_2)
      
      matches_df <- chunk_df %>%
        dplyr::mutate(index_number = 1:dplyr::n()) 
      
      n_total_reads <- n_total_reads + length(forward_reads_string_set)
      
      if(nrow(matches_df) == 0) {
        next
      } else if (lp > length(cumulative_count_df)) { 
        cumulative_count_df[[lp]] = (matches_df %>% 
                                     dplyr::count(index_1, index_2, forward_read_cl_barcode))
      } else {
        cumulative_count_df[[lp]] = (matches_df %>% 
          dplyr::count(index_1, index_2, forward_read_cl_barcode) %>% rbind(cumulative_count_df[[lp]])) # new
      }
    }
    close(forward_stream)
    
    lp = lp + 1 # new
    
    close(index_2_stream)
    close(index_1_stream)
    
    if (!is.na(write_interval) & (i %% write_interval == 0)) {
      print(paste0('saving cumulative count df at iteration, ', i))
      out = cumulative_count_df %>% dplyr::bind_rows() # new
      #saveRDS(cumulative_count_df, 'temporary_cumulative_count_df.Rds')
      saveRDS(out,  paste(save_loc,'temporary_cumulative_count_df.Rds', sep='/')) # new
    }
  }
  
  print('saving final cumulative_count_df ')
  out = cumulative_count_df %>% dplyr::bind_rows() # new
  #saveRDS(cumulative_count_df, 'temporary_cumulative_count_df.Rds')
  saveRDS(out, paste(save_loc, 'cumulative_count_df.Rds', sep = '/')) # new
  
  cumulative_count_df <- cumulative_count_df %>%
    dplyr::bind_rows() %>% # new 
    dplyr::group_by(index_1, index_2, forward_read_cl_barcode) %>%
    dplyr::summarise(n = sum(n, na.rm = T)) %>%
    dplyr::ungroup()
  
  return(cumulative_count_df)
}

# this function is for DRAGEN-formatted files
write_df_from_fastq_DRAGEN <- function(
    forward_read_fastq_files,
    CL_BC_LENGTH = 24,
    PLATE_BC_LENGTH = 8,
    WELL_BC_LENGTH = 8,
    save_loc = NULL
) 
  {
  require(tidyverse)
  require(magrittr)

  cumulative_count_df_uncollapsed <- list()
  
  lp = 1
  for (i in 1:length(forward_read_fastq_files)) {
    forward_stream <- ShortRead::FastqStreamer(forward_read_fastq_files[i])
    
    print(paste0("forward read file is: ", forward_read_fastq_files[i]))
    
    file_name_without_path <- word(forward_read_fastq_files[i],-1,sep=fixed("/")) ## extract file name without directoryname
    flow_cell <- word(file_name_without_path,1,sep=fixed("_")) ## get flow_cell 
    flow_lane <- word(file_name_without_path,2,sep=fixed("_")) ## get flow_lane
    
    j <- 0
    repeat{
      forward_reads_chunk <- ShortRead::yield(forward_stream)
      
      if (length(forward_reads_chunk) == 0) {break}
      
      print(paste0('Processing read ', j * 1e6, ' to ', 
                   j*1e6 + length(forward_reads_chunk), 
                   ' from file ', i)) #(ShortRead)The default size for both streams and samples is 1M records;
      j <- j + 1
      
      forward_reads_string_set <- forward_reads_chunk %>%
        ShortRead::sread() %>%
        Biostrings::DNAStringSet()
      
      forward_reads_cl_barcode <- XVector::subseq(forward_reads_string_set, 1, CL_BC_LENGTH) %>%
        as.character()
      
      read_header_str <- as.character(forward_reads_chunk@id)
      read_header_len <- nchar(read_header_str)
      
      print("Return the counts for each barcode and index pair in the chunk")
      cumulative_count_df_uncollapsed[[lp]] <- data.frame(indeces = substr(x = read_header_str, 
                                                               (read_header_len)-(PLATE_BC_LENGTH+WELL_BC_LENGTH),
                                                               read_header_len),
                                              forward_read_cl_barcode = forward_reads_cl_barcode,
                                              flowcell_name=flow_cell,
                                              flowcell_lane=flow_lane)  %>%
        dplyr::count(indeces, forward_read_cl_barcode, flowcell_name, flowcell_lane) 
      
      lp <- lp + 1
    }
    close(forward_stream)
  }
  
  print('saving final cumulative_count_df ')
  ## cumulative counts separate across different flowcells and flowcell lanes.
  cumulative_count_df_uncollapsed = cumulative_count_df_uncollapsed %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(indeces, forward_read_cl_barcode, flowcell_name, flowcell_lane) %>% 
    dplyr::summarise(n = sum(n, na.rm = T)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(index_1 = word(indeces, 1, sep = fixed("+")),
                  index_2 = word(indeces, 2, sep = fixed("+"))) %>% 
    dplyr::select(-indeces)
  
  ## get counts summed across different flow cells and lanes. Summing across lanes was implicit before and implicit in other sequencer formats
  cumulative_count_collapsed_across_flowcells_df <- cumulative_count_df_uncollapsed %>% 
    dplyr::group_by(index_1, index_2, forward_read_cl_barcode) %>%
    dplyr::summarise(n = sum(n, na.rm = T)) %>% 
    dplyr::ungroup()
  print (paste("collapsed across", length(cumulative_count_df_uncollapsed$flowcell_name %>% unique()), "flowcells"))
  if(!is.null(save_loc)){
    write_csv(cumulative_count_df_uncollapsed, file =  paste0(save_loc, '/raw_counts_uncollapsed.csv')) 

  }
  return(cumulative_count_collapsed_across_flowcells_df) ## we want the collapsed data  
}
