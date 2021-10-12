suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))

#source("../src/load_libraries.R")
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr)) #write_delim 
suppressPackageStartupMessages(library(stringr)) #str_detect
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(tidyr)) #pivot_wider
suppressPackageStartupMessages(library(reshape2))

# source("../src/write_df_from_fastq.R", local=T)

## write_df_from_fastq
## takes PRISM miseq or hiseq fastqs and returns a dataframe
## with the number of barcodes from each cell line in each well
##
## takes: 
##     forward_read_fastq_files - a vector of fastq file paths.
##     index_1_files - a vector of fastq file paths.
##     index_2_files - a vector of fastq file paths. 
##     write_interval - integer specifying how often a temporary count file should be written. 
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
        dplyr::mutate(index_number = 1:n()) 
      
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

## print_args
## writes configuration to file
##
## takes: 
##      args: args object from argparse
print_args <- function(args){
  config <- data.frame(args=names(args), values=unname(unlist(args)))
  config_path = paste(
    args$out, 
    "config.txt",
    sep="/"
  )
  print(paste("Saving config.txt file in :", config_path))
  write_delim(config, config_path, delim = ": ", col_names=F)
}

# create parser object
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")

parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-f", "--fastq", default="fastq/",
                    help="Path to directory containing fastq files")
parser$add_argument("-i1", "--index_1", default="", help = "Index 1 code")
parser$add_argument("-i2", "--index_2", default="", help = "Index 2 code")
parser$add_argument("-b", "--barcode_suffix", default="", help = "Barcode Read Files code")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

#Save arguments
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

read_directory_contents <- c(args$fastq) %>% 
  purrr::map(list.files, full.names = T) %>%
  purrr::reduce(union)

barcode_read_files <- read_directory_contents %>%
  purrr::keep(str_detect, fixed(args$barcode_suffix)) %>%
  sort()

# plates
index_1_files <- read_directory_contents %>%
  purrr::keep(str_detect, fixed(args$index_1)) %>%
  sort()

# wells
index_2_files <- read_directory_contents %>%
  purrr::keep(str_detect, fixed(args$index_2)) %>%
  sort()

print(paste("num index_1 files", length(index_1_files)))
print(paste("num index_2 files", length(index_2_files)))
print(paste("num barcode files", length(barcode_read_files)))

print("creating read count file")
raw_counts <- write_df_from_fastq(forward_read_fastq_files = barcode_read_files,
                                  index_1_file = index_1_files,
                                  index_2_file = index_2_files,
                                  write_interval = 100)


rc_out_file = paste(
  args$out,
  'raw_counts.csv',
  sep='/'
)
print(paste("writing to file: ", rc_out_file))
write.csv(raw_counts, rc_out_file, row.names=F, quote=F)

## Make GCTX file
## For compatibility with other tools

# 
# counts_df = with(
#   filtered_counts[!is.na(filtered_counts$LUA),], #Remove rows where cell line LUA is not known
#   data.frame(
#     cid = profile_id,
#     rid=LUA,
#     value=n )
# )   #Extract matrix columns
# 
# #Pivot 
# counts_mat <- counts_df %>% 
#   pivot_wider(
#     names_from = cid,
#     id_cols = rid,
#     values_from = value,
#     values_fn=as.numeric
#   )
# 
# rids = counts_mat$rid
# counts_mat <- as.matrix(counts_mat[-1])
# rownames(counts_mat) <- rids
# 
# rownames(sample_meta) <- sample_meta$profile_id
# rownames(cell_line_meta) <- cell_line_meta$LUA
# 
# col_desc = sample_meta[sample_meta$profile_id %in% colnames(counts_mat),]
# row_desc = cell_line_meta[cell_line_meta$LUA %in% rids,]
# 
# #counts_gct <- new("GCT", mat=counts_mat, rid = rownames(counts_mat), cid = colnames(counts_mat), rdesc=row_desc, cdesc=col_desc)
# 
# filtrc_out_gct = paste(
#   args$out,
#   'Level2_counts',
#   sep='/'
# )
# #print(paste0("Writing GCTx file to ", filtrc_out_gct))
# #write_gctx(counts_gct, filtrc_out_gct)

