options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
source("./src/compute_l2fc.R")
source("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("-ct", "--control_type", default="negcon", help="pert_type to use as control")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--ctrl_cols", default="cell_set,day", 
                    help = "columns used to collapse controls to generate l2fc")
parser$add_argument("--cell_line_cols", default="project_code,depmap_id,ccle_name", 
                    help = "Columns that can describe a cell line")
parser$add_argument("-ccn", "--count_col_name", default="normalized_n", 
                    help = "column containing counts with which to calculate l2fc")
parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
parser$add_argument("-o","--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up parameters and run compute_l2fc ----
control_type = args$control_type
normalized_counts= data.table::fread(args$normalized_counts, header= TRUE, sep= ',', data.table= FALSE)
sig_cols = unlist(strsplit(args$sig_cols, ","))
ctrl_cols = unlist(strsplit(args$ctrl_cols, ","))
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))
count_col_name = args$count_col_name
count_threshold = as.numeric(args$count_threshold)

print("Collapsing tech reps and computing log-fold change ...")
l2fc= compute_l2fc(normalized_counts= normalized_counts, 
                   control_type= control_type, 
                   sig_cols= sig_cols, 
                   ctrl_cols= ctrl_cols, 
                   count_col_name= count_col_name, 
                   count_threshold = count_threshold,
                   cell_line_cols= cell_line_cols)

# Write out file ----
l2fc_outpath= paste(args$out, "l2fc.csv", sep= "/")
print(paste0('Writing out l2fc file to ', l2fc_outpath))
write.csv(l2fc, l2fc_outpath, row.names= FALSE, quote= FALSE)

# Ensure that l2fc file was successfully generated ----
check_file_exists(l2fc_outpath)
