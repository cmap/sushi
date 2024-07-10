suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(magrittr))
source("/workspace/R/collapse_counts.R")
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--lfc", default="l2fc.csv",
                    help="path to file containing l2fc values")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

lfc_values = read.csv(args$lfc)

print("collapsing counts")
collapsed_counts = collapse_counts(lfc_values)

collapsed_count_out_file = paste(
  args$out,
  "collapsed_values.csv",
  sep='/'
)

write.csv(collapsed_counts, collapsed_count_out_file, row.names=F, quote=F)
