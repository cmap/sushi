options(cli.unicode = FALSE)
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
suppressPackageStartupMessages(library(WGCNA))
source("/workspace/scripts/src/QC_images.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("-c", "--uncollapsed_raw_counts", default="raw_counts_uncollapsed.csv",
                    help="path to file containing uncollapsed raw counts file")
parser$add_argument("--raw_counts", default= "raw_counts.csv", help="path to raw counts file")
parser$add_argument("--annotated_counts", default= "annotated_counts.csv",
                    help= "path to file containing annotated counts")
parser$add_argument("--filtered_counts", default= "filtered_counts.csv", help= "path to filtered_counts file")
parser$add_argument("--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("--l2fc", default="l2fc.csv", help= "path to l2fc file")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help = "Sample metadata")
parser$add_argument("--CB_meta", default="/data/CB_meta.csv", help = "control barcode metadata")
parser$add_argument("--cell_set_meta", default="cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("--cell_line_cols", default= 'DepMap_ID,CCLE_name',
                    help= "Columns that identify cell lines or barcodes")
parser$add_argument("--id_cols", default= 'pcr_plate,pcr_well',
                    help= "Columns to identify each PCR well")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help= 'Columns used to identify the treatment conditions')
parser$add_argument("--control_type", default = "negcon",
                    help= "how negative control wells are distinguished in the trt_type column")
parser$add_argument("--count_threshold", default=40, help= "Low counts threshold")
parser$add_argument("--reverse_index2", default=FALSE, help = "Reverse index 2")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}

# Read in files and pull out parameters ----
# Pipeline outputs
raw_counts_uncollapsed= data.table::fread(args$uncollapsed_raw_counts, header= TRUE, sep= ',')
raw_counts= data.table::fread(args$raw_counts, header= TRUE, sep= ',')
annotated_counts= data.table::fread(args$annotated_counts, header= TRUE, sep= ',')
filtered_counts= data.table::fread(args$filtered_counts, header= TRUE, sep= ',')
if(file.exists(args$normalized_counts)) {
  normalized_counts= data.table::fread(args$normalized_counts, header=TRUE, sep=',', data.table=FALSE)
} else {
  normalized_counts= NA
}
l2fc= data.table::fread(args$l2fc, header= TRUE, sep= ',')

# Metadata files
sample_meta= data.table::fread(args$sample_meta, header= TRUE, sep= ',', data.table= FALSE)
CB_meta= data.table::fread(args$CB_meta, header=TRUE, sep=',', data.table=FALSE)
cell_set_meta = data.table::fread(args$cell_set_meta, header=TRUE, sep=',', data.table=FALSE)

# Parameters
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
id_cols= unlist(strsplit(args$id_cols, ","))
sig_cols= unlist(strsplit(args$sig_cols, ","))
control_type = args$control_type
count_threshold= as.numeric(args$count_threshold)

# # If flag passed, use cell_set_meta file generated for the project via CellDB
# if (args$db_flag) {
#   print("Calling cell_set_meta generated using CellDB")
#   cell_set_meta = read.csv("cell_set_meta.csv")
#   # Otherwise, use static file
# } else {
#   print("Using static cell set metadata file to locate cell information.")
#   cell_set_meta = read.csv(args$cell_set_meta)
# }

print("Generating QC images ...")
QC_images(raw_counts_uncollapsed= raw_counts_uncollapsed, 
          raw_counts= raw_counts, 
          annotated_counts= annotated_counts, 
          filtered_counts= filtered_counts, 
          normalized_counts= normalized_counts, 
          l2fc= l2fc, 
          sample_meta= sample_meta, 
          CB_meta= CB_meta, 
          cell_set_meta= cell_set_meta,
          cell_line_cols= c('DepMap_ID', 'CCLE_name'), 
          id_cols= id_cols, 
          sig_cols= sig_cols,
          control_type= control_type, 
          count_threshold= count_threshold, 
          reverse_index2= args$reverse_index2,
          out= args$out)
