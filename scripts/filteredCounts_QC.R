suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(grDevices))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(scales)) # for out of bound handling in plots
suppressPackageStartupMessages(library(ggpmisc)) # with ggplot to add fit line and labels
source("/workspace/R/QC_images.R")

# Parser ----
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("--raw_counts", default= "raw_counts.csv", help="path to file containing raw counts")
parser$add_argument("--annotated_counts", default="annotated_counts.csv",
                    help="path to file containing annotated counts")
# parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
#                     help="path to file containing filtered counts")
parser$add_argument("--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("--CB_meta", default="/data/CB_meta.csv", help = "control barcode metadata")
parser$add_argument("--cell_set_meta", default="cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help = "columns used to generate signature ids")
parser$add_argument("--count_col_name", default="normalized_n", 
                    help = "column containing counts with which to calculate l2fc")
parser$add_argument("--count_threshold", default=40, 
                    help = "Low counts threshold")
parser$add_argument("--reverse_index2", default=FALSE, help = "Reverse index 2")
parser$add_argument("--control_type", default = "negcon",
                    help = "how negative control wells are distinguished in the trt_type column")
# parser$add_argument("--db_flag", action="store_true", default=FALSE, help = "Use CellDB to locate cell set information")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

# Read in files and pull out parameters ----
sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',', data.table= F)
raw_counts= data.table::fread(args$raw_counts, header= T, sep= ',', data.table= F)
annotated_counts= data.table::fread(args$annotated_counts, header= T, sep= ',', data.table= F)
if(file.exists(args$normalized_counts)) {
  normalized_counts= data.table::fread(args$normalized_counts, header= T, sep= ',', data.table= F)
} else {
  normalized_counts= NA
}
CB_meta= read.csv(args$CB_meta)
sig_cols = unlist(strsplit(args$sig_cols, ","))
count_col_name = args$count_col_name
count_threshold_arg= args$count_threshold
count_threshold = as.numeric(count_threshold_arg)
cell_set_meta = read.csv(args$cell_set_meta)
control_type = args$control_type

# # If flag passed, use cell_set_meta file generated for the project via CellDB
# if (args$db_flag) {
#   print("Calling cell_set_meta generated using CellDB")
#   cell_set_meta = read.csv("cell_set_meta.csv")
#   # Otherwise, use static file
# } else {
#   print("Using static cell set metadata file to locate cell information.")
#   cell_set_meta = read.csv(args$cell_set_meta)
# }

print("generating filtered counts QC images")
#QC_images(annotated_counts, filtered_counts, normalized_counts,
#          CB_meta, cell_set_meta, args$out, sig_cols, count_col_name)
QC_images(sample_meta= sample_meta,
          annotated_counts= annotated_counts, 
          raw_counts= raw_counts,
          normalized_counts= normalized_counts,
          CB_meta= CB_meta, 
          cell_set_meta= cell_set_meta, 
          control_type = control_type,
          out= args$out, 
          sig_cols= sig_cols, 
          count_col_name= count_col_name, 
          count_threshold= count_threshold,
          reverse_index2= args$reverse_index2)
