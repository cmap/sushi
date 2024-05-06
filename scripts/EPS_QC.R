# Generates the QC table for EPS 

##Assay pool information will need to be added when merging has been implemented in pipeline
## so now QC table has cell set information rather than assay pool information
# import necessary libraries and functions
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(data.table))

#---- Read arguments ----
# initialize parser
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-b", "--base_dir", default="", help="Input Directory with data from whole screen")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")
parser$add_argument("-n", "--name", default="", help = "Build name. Default is none")
# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

base_dir <- args$base_dir
out_dir <- args$out
build_name <- args$name

count_threshold <- 40 ## could be passed as an option in the future


if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = T)} ## create output dir if it doesn't exist

## read data
norm_count_files <- list.files(base_dir, pattern = "normalized_counts.csv", full.names = T)
if (length(norm_count_files) == 1) {
    norm_counts <- data.table::fread(norm_count_files[[1]])
} else {
    stop(paste("There are", length(norm_count_files), "normalized count files in this directory. Please ensure there is one and try again."),
         call. = FALSE)
}


## get negcon counts at day 10
cell_counts_negcon <- norm_counts %>% 
    dplyr::filter(!is.na(CCLE_name)) %>% 
    dplyr::filter(day==10) %>%
    dplyr::filter(trt_type=="negcon")


## get median counts in negcon and apply count threshold to generate the QC table
med_cell_counts_qc <- cell_counts_negcon %>% 
    dplyr::group_by(CCLE_name, DepMap_ID, cell_set, day, pert_plate) %>% 
    dplyr::summarize(med_log_raw_counts=median(log2_n), 
                     med_log_norm_counts= median(log2_normalized_n),
                     med_raw_counts= median(n)) %>% 
    dplyr::mutate(pass_raw_count_qc=med_raw_counts>count_threshold) %>% 
    dplyr::ungroup() 

## in log2 space, the threshold should be at log2(40+threshold)

## save the output along with build name? or without ?
write_csv(med_cell_counts_qc, paste0(out_dir, "/", build_name, "EPS_QC_TABLE.csv"))











