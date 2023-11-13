suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(grDevices))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(prismSeqR))
suppressPackageStartupMessages(library(scales)) # for out of bound handling in plots
suppressPackageStartupMessages(library(ggpmisc)) # with ggplot to add fit line and labels

parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("--annotated_counts", default="annotated_counts.csv",
                    help="path to file containing annotated counts")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "control barcode metadata")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help = "columns used to generate signature ids")
parser$add_argument("--count_col_name", default="normalized_n", 
                    help = "column containing counts with which to calculate l2fc")
parser$add_argument("--count_threshold", default=40, 
                    help = "Low counts threshold")
parser$add_argument("--db_flag", action="store_true", default=TRUE, help = "Use CellDB to locate cell set information")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

sample_meta = read.csv(args$sample_meta)
annotated_counts= read.csv(args$annotated_counts)
filtered_counts = read.csv(args$filtered_counts)
if(file.exists(args$normalized_counts)) {
  normalized_counts= read.csv(args$normalized_counts)
} else {
  normalized_counts= NA
}
CB_meta= read.csv(args$CB_meta)
sig_cols = unlist(strsplit(args$sig_cols, ","))
count_col_name = args$count_col_name
count_threshold_arg= args$count_threshold
count_threshold = as.numeric(count_threshold_arg)

# If flag passed, use cell_set_meta file generated for the project via CellDB
if (args$db_flag) {
  print("Calling cell_set_meta generated using CellDB")
  cell_set_meta = read.csv("cell_set_meta.csv")
  # Otherwise, use static file
} else {
  print("Using static cell set metadata file to locate cell information.")
  cell_set_meta = read.csv("../metadata/cell_set_meta.csv")
}

print("generating filtered counts QC images")
#QC_images(annotated_counts, filtered_counts, normalized_counts,
#          CB_meta, cell_set_meta, args$out, sig_cols, count_col_name)
QC_images(sample_meta= sample_meta,
          annotated_counts= annotated_counts, 
          filtered_counts= filtered_counts,
          normalized_counts= normalized_counts,
          CB_meta= CB_meta, cell_set_meta= cell_set_meta, 
          out= args$out, 
          sig_cols= sig_cols, count_col_name= count_col_name, count_threshold= count_threshold)
