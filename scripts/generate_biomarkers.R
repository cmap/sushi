library(argparse)
library(dplyr)
library(scam)
library(magrittr)
library(tidyr)
library(reshape2)
library(tibble)
library(cdsrbiomarker)
source("/workspace/R/generate_biomarkers.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--collapsed_values", default="collapsed_values.csv",
                    help="path to file containing collapsed l2fc values")
parser$add_argument("--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == ""){
  args$out = args$wkdir
}

collapsed_values = read.csv(args$collapsed_values)

print("generating biomarker tables")
biomarker_out = generate_biomarkers(collapsed_values)

lin_table = biomarker_out$lin_table
rf_table = biomarker_out$rf_table
disc_table = biomarker_out$disc_table

lin_out = paste(args$out, "lin_table.csv", sep='/')
rf_out = paste(args$out, "rf_table.csv", sep='/')
disc_out = paste(args$out, "disc_table.csv", sep='/')

print("writing out biomarker tables")
write.csv(lin_table, lin_out, row.names=F, quote=F)
write.csv(rf_table, rf_out, row.names=F, quote=F)
write.csv(disc_table, disc_out, row.names=F, quote=F)
