suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
#suppressPackageStartupMessages(library(scam))
#suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(reshape2))
library(prismSeqR)

#Need Arguments
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("-ct", "--control_type", default="negcon", help="trt_type to use as control")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help = "columns used to generate signature ids")
parser$add_argument("-ccn", "--count_col_name", default="normalized_n", 
                    help = "column containing counts with which to calculate l2fc")
parser$add_argument("-o","--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

control_type = args$control_type
normalized_counts = read.csv(args$normalized_counts)
sig_cols = unlist(strsplit(args$sig_cols, ","))
count_col_name = args$count_col_name

print("computing log-fold change")
l2fc = compute_l2fc(normalized_counts, control_type, sig_cols, count_col_name)

l2fc_out = paste(args$out, "l2fc.csv", sep="/")
write.csv(l2fc, l2fc_out, row.names=F, quote=F)
