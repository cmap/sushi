library(argparse)
library(magrittr)
library(tidyverse)
library(data.table)
source ("./src/qc_functions.R")

# Argument parser ----
parser <- ArgumentParser()
parser$add_argument("--raw_counts_uncollapsed", default= "raw_counts_uncollapsed.csv",
                    help= "path to file containing uncollapsed raw counts file")
parser$add_argument("--sample_meta", default="sample_meta.csv", help= "Sample metadata")
parser$add_argument("--cell_line_meta", default="cell_line_meta.csv", help= "Cell line metadata")
parser$add_argument("--CB_meta", default= "CB_meta.csv", help= "Control Barcode metadata")
parser$add_argument("--l2fc", default= "l2fc.csv", help= "l2fc file")
parser$add_argument("--collapsed_l2fc", default= "collapsed_l2fc.csv", help= "collapsed l2fc file")
parser$add_argument("--normalized_counts", default= "normalized_counts.csv", help= "normalized counts file")