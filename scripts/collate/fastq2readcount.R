options(cli.unicode = FALSE)
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse)) ###### debug
suppressPackageStartupMessages(library(magrittr))
source("write_df_from_fastq.R")

# Argument parser ----
parser <- argparse::ArgumentParser()
# specify desired options
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
parser$add_argument("-w", "--write_interval", default=500, help = "integer for how often a temp count file is written.")
parser$add_argument("-l", "--barcode_lengths", default='8,8,24',
                    help = "Three number, comma-separated string that denotes, PLATE_BC_LEN,WELL_BC_LEN,CELL_BC_LEN, e.g. 8,8,24")
parser$add_argument("-s", "--seq_type", default = "", help = "Designate DRAGEN if from DRAGEN")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

#if args$out doesn't exist, create it
if (!dir.exists(args$out)){
  dir.create(args$out)
}

read_directory_contents <- c(args$fastq) %>% 
  purrr::map(list.files, full.names = T) %>%
  purrr::reduce(union)

barcode_read_files <- read_directory_contents %>%
  purrr::keep(stringr::str_detect, fixed(args$barcode_suffix)) %>%
  sort()

if(args$seq_type == "DRAGEN"){
  barcode_read_files %<>% purrr::keep(stringr::str_detect, fixed("_R1_001.fastq"))
} else {
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
}

print(paste("num barcode files", length(barcode_read_files)))

print("creating read count file")

if (is.null(args$write_interval)){
  write_interval = NA
} else {
  write_interval = as.numeric(args$write_interval)
}

#barcode_read_lengths
read_lengths = as.integer(unlist(strsplit(args$barcode_lengths, ",")))

if (args$seq_type == "DRAGEN"){
  raw_counts <- write_df_from_fastq_DRAGEN(forward_read_fastq_files = barcode_read_files,
                                       PLATE_BC_LENGTH = read_lengths[1],
                                       WELL_BC_LENGTH = read_lengths[2],
                                       CL_BC_LENGTH = read_lengths[3],
                                       save_loc =  args$out)
}else{
  raw_counts <- write_df_from_fastq(forward_read_fastq_files = barcode_read_files,
                                    index_1_file = index_1_files,
                                    index_2_file = index_2_files,
                                    write_interval = write_interval, 
                                    PLATE_BC_LENGTH = read_lengths[1],
                                    WELL_BC_LENGTH = read_lengths[2],
                                    CL_BC_LENGTH = read_lengths[3],
                                    save_loc = args$out)
}




rc_out_file = paste(
  args$out,
  'raw_counts.csv',
  sep='/'
)
print(paste("writing to file: ", rc_out_file))
write.csv(raw_counts, rc_out_file, row.names=F, quote=F)

