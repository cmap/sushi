options(cli.unicode = FALSE)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(grDevices))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales)) # for out of bound handling in plots
suppressPackageStartupMessages(library(ggpmisc)) # with ggplot to add linear fit labels
suppressPackageStartupMessages(library(WGCNA)) # for faster correlations
source("/workspace/scripts/src/QC_images.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE, help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("-c", "--raw_counts_uncollapsed", default="raw_counts_uncollapsed.csv",
                    help="path to file containing uncollapsed raw counts file")
parser$add_argument("--raw_counts", default= "raw_counts.csv", help="path to raw counts file")
parser$add_argument("--annotated_counts", default= "annotated_counts.csv",
                    help= "path to file containing annotated counts")
parser$add_argument("--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("--lfc", default="l2fc.csv", help= "path to l2fc file")
parser$add_argument("--cell_line_cols", default= 'DepMap_ID,CCLE_name',
                    help= "Columns that identify cell lines or barcodes")
parser$add_argument("--id_cols", default= 'pcr_plate,pcr_well',
                    help= "Columns to identify each PCR well")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help= 'Columns used to identify the treatment conditions')
parser$add_argument("--control_type", default = "negcon",
                    help= "how negative control wells are distinguished in the trt_type column")
parser$add_argument("--count_threshold", default=40, help= "Low counts threshold")
parser$add_argument("--reverse_index2", type="logical", default=FALSE,
                    help= "Reverse complement of index 2 for NovaSeq and NextSeq")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}

# Read in files and pull out parameters ----
sample_meta= data.table::fread(args$sample_meta, header= TRUE, sep= ',')

# Pipeline outputs
raw_counts_uncollapsed= data.table::fread(args$raw_counts_uncollapsed, header= TRUE, sep= ',')
raw_counts= data.table::fread(args$raw_counts, header= TRUE, sep= ',')
annotated_counts= data.table::fread(args$annotated_counts, header= TRUE, sep= ',')
if(file.exists(args$normalized_counts)) {
  normalized_counts= data.table::fread(args$normalized_counts, header=TRUE, sep=',', data.table=FALSE)
} else {
  normalized_counts= NA
}
l2fc= data.table::fread(args$lfc, header= TRUE, sep= ',')

# Parameters
cell_line_cols = unlist(strsplit(args$cell_line_cols, ","))
id_cols= unlist(strsplit(args$id_cols, ","))
sig_cols= unlist(strsplit(args$sig_cols, ","))
control_type = args$control_type
count_threshold= as.numeric(args$count_threshold)
#

# Call QC images function ----
print("Calling QC images ...")
QC_images(raw_counts_uncollapsed= raw_counts_uncollapsed, 
          raw_counts= raw_counts, 
          annotated_counts= annotated_counts, 
          normalized_counts= normalized_counts, 
          l2fc= l2fc, 
          sample_meta= sample_meta,
          cell_line_cols= c('DepMap_ID', 'CCLE_name'), 
          id_cols= id_cols, 
          sig_cols= sig_cols,
          control_type= control_type, 
          count_threshold= count_threshold, 
          reverse_index2= args$reverse_index2,
          out= args$out)
