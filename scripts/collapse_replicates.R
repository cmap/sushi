suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(magrittr))
library(prismSeqR)

parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", dest="verbose", 
                    help="Print little output")
parser$add_argument("-c", "--lfc", default="l2fc.csv",
                    help="path to file containing l2fc values")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help = "columns used to identify a unique condition")
parser$add_argument("--cell_line_cols", default="project_code,DepMap_ID,CCLE_name", 
                    help = "Columns that can describe a cell line")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

lfc_values= data.table::fread(args$lfc, header=T, sep=',')

print("Collapsing biological replicates ...")
sig_cols= unlist(strsplit(args$sig_cols, ","))
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))
collapsed_l2fc= collapse_bio_reps(lfc_values, sig_cols, cell_line_cols)

# Write out file
collapsed_l2fc_path= paste(args$out, "collapsed_l2fc.csv", sep='/')
write.csv(collapsed_l2fc, collapsed_l2fc_path, row.names=F, quote=F)