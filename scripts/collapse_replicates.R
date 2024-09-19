options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(magrittr))
source("./src/collapse_bio_reps.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action= "store_true", default= TRUE,
                    help= "Print extra output [default]")
parser$add_argument("-q", "--quietly", action= "store_false", dest= "verbose", 
                    help= "Print little output")
parser$add_argument("-c", "--lfc", default= "l2fc.csv",
                    help= "path to file containing l2fc values")
parser$add_argument("--sig_cols", default= "cell_set,treatment,dose,dose_unit,day", 
                    help= "columns used to identify a unique condition")
parser$add_argument("--cell_line_cols", default= "project_code,depmap_id,ccle_name", 
                    help= "Columns that can describe a cell line")
parser$add_argument("-o", "--out", default= getwd(), help= "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Collapse biological replicates ----
lfc_values= data.table::fread(args$lfc, header=T, sep=',', data.table=F)
sig_cols= unlist(strsplit(args$sig_cols, ","))
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))

print("Collapsing biological replicates ...")
collapsed_l2fc= collapse_bio_reps(l2fc= lfc_values, sig_cols= sig_cols, cell_line_cols= cell_line_cols)

# Write out file ----
collapsed_l2fc_outpath= paste(args$out, 'collapsed_l2fc.csv', sep='/')
print(paste0('Writing out collapsed l2fc file to ', collapsed_l2fc_outpath))
write.csv(x= collapsed_l2fc, file= collapsed_l2fc_outpath, row.names= FALSE, quote= FALSE)
