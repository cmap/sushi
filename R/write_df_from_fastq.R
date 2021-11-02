#' write_df_from_fastq
#' 
#' takes PRISM miseq or hiseq fastqs and returns a dataframe
#' with the number of barcodes from each cell line in each well
#' 
#' 
#' 
#' @param forward_read_fastq_files Vector of fastq file paths
#' @param index_1_files Vector of fastq file paths
#' @param index_2_files Vector of fastq file paths
#' @param write_interval integer for how often a temp count file is written. 
#' @return cumulative_count_df A data.frame of readcounts by index_1, index_2 and forward_read_cl_barcode
#' @export 
write_df_from_fastq <- function(
  forward_read_fastq_files, 
  index_1_files, 
  index_2_files,
  write_interval = NA) 
{
  
  CL_BC_LENGTH <- 24
  PLATE_BC_LENGTH <- 8
  WELL_BC_LENGTH <- 8
  FORWARD_CONSTANT_REGION_1_LENGTH <- 34
  FORWARD_CONSTANT_REGION_2_LENGTH <- 7
  
  n_total_reads <- 0
  cumulative_count_df <- data.frame()
  
  for (i in 1:length(forward_read_fastq_files)) {
    forward_stream <- ShortRead::FastqStreamer(forward_read_fastq_files[i])
    index_1_stream <- ShortRead::FastqStreamer(index_1_files[i])
    index_2_stream <- ShortRead::FastqStreamer(index_2_files[i])
    
    j <- 0
    repeat{
      forward_reads_chunk <- ShortRead::yield(forward_stream)
      index_1_chunk <- ShortRead::yield(index_1_stream)
      index_2_chunk <- ShortRead::yield(index_2_stream)
      
      if (length(forward_reads_chunk) == 0) {break}
      
      print(paste0('Processing read ', j * 1e6, ' to ', 
                   j*1e6 + length(forward_reads_chunk), 
                   ' from file ', i))
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
      
      index_1 <- XVector::subseq(index_1_string_set, 1, WELL_BC_LENGTH)
      index_2 <- XVector::subseq(index_2_string_set, 1, PLATE_BC_LENGTH)
      
      forward_reads_cl_barcode = as.character(forward_reads_cl_barcode)
      index_1 = as.character(index_1)
      index_2 = as.character(index_2)
      chunk_df <- data.frame(forward_read_cl_barcode = forward_reads_cl_barcode, index_1 = index_1, index_2 = index_2)
      
      matches_df <- chunk_df %>%
        dplyr::mutate(index_number = 1:dplyr::n()) 
      
      n_total_reads <- n_total_reads + length(forward_reads_string_set)
      
      if(nrow(matches_df) == 0) {
        next
      } else {
        cumulative_count_df <- matches_df %>%
          dplyr::count(index_1, index_2, forward_read_cl_barcode) %>%
          rbind(cumulative_count_df)
      }
    }
    close(forward_stream)
    
    close(index_2_stream)
    close(index_1_stream)
    
    if (!is.na(write_interval) & (i %% write_interval == 0)) {
      print(paste0('saving cumulative count df at iteration, ', i))
      saveRDS(cumulative_count_df, 'temporary_cumulative_count_df.Rds')
    }
  }
  
  cumulative_count_df <- cumulative_count_df %>%
    dplyr::group_by(index_1, index_2, forward_read_cl_barcode) %>%
    dplyr::summarise(n = sum(n, na.rm = T)) %>%
    dplyr::ungroup()
  
  return(cumulative_count_df)
}
