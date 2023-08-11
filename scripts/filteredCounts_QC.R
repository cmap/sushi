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

parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("--annotated_counts", default="annotated_counts.csv",
                    help="path to file containing annotated counts")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "control barcode metadata")
parser$add_argument("--cell_set_meta", default="../metadata/cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help = "columns used to generate signature ids")
parser$add_argument("-ccn", "--count_col_name", default="normalized_n", 
                    help = "column containing counts with which to calculate l2fc")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

annotated_counts= read.csv(args$annotated_counts)
filtered_counts = read.csv(args$filtered_counts)
normalized_counts= read.csv(args$normalized_counts)
CB_meta= ead.csv(args$CB_meta)
cell_set_meta = read.csv(args$cell_set_meta)
sig_cols = unlist(strsplit(args$sig_cols, ","))
count_col_name = args$count_col_name


print("generating filtered counts QC images")
QC_images(annotated_counts, filtered_counts, normalized_counts,
          CB_meta, cell_set_meta, args$out, sig_cols, count_col_name)
