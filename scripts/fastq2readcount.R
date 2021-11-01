source("./R/write_df_from_fastq.R")

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
parser <- argparse::ArgumentParser()
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
  purrr::keep(stringr::str_detect, fixed(args$barcode_suffix)) %>%
  sort()

# plates
index_1_files <- read_directory_contents %>%
  purrr::keep(stringr::str_detect, fixed(args$index_1)) %>%
  sort()

# wells
index_2_files <- read_directory_contents %>%
  purrr::keep(stringr::str_detect, fixed(args$index_2)) %>%
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

