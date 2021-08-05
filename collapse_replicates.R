suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tibble))

## collapse_counts
## collapses filtered normalized counts and computes MAD/sqrt(n) metrics.
## cell lines with MAD/sqrt(n) > 0.5 are filtered out, and file with filtered out cell lines is written. 
## log10(median counts) vs. MAD/sqrt(n) graph is saved, and collapsed filtered count table is returned
##
## takes:
##      filtered_normalized_counts - normalized counts with bad replicates already filtered out
collapse_counts = function(l2fc) {
  collapsed_counts = l2fc %>% 
    filter(control_pass_QC) %>% 
    group_by_at(setdiff(names(.), c("bio_rep", "sum_normalized_n", "control_mad_sqrtN", "l2fc", "control_pass_QC", "control_median_normalized_n"))) %>% 
    dplyr::summarise(trt_median_normalized_n = median(sum_normalized_n),
                     median_l2fc = median(l2fc),
                     trt_mad_sqrtN = mad(log10(sum_normalized_n))/sqrt(n())) %>% 
    ungroup() %>% 
    mutate(trt_pass_QC = ifelse(trt_mad_sqrtN > 0.5, F, T)) %>% 
    dplyr::relocate(trt_median_normalized_n, trt_mad_sqrtN, trt_pass_QC, median_l2fc, .after=last_col())
  
  return(collapsed_counts)
}


parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--lfc", default="l2fc.csv",
                    help="path to file containing l2fc values")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

lfc_values = read.csv(args$lfc)

print("collapsing counts")
collapsed_counts = collapse_counts(lfc_values)

collapsed_count_out_file = paste(
  args$out,
  "collapsed_values.csv",
  sep='/'
)

write.csv(collapsed_counts, collapsed_count_out_file, row.names=F, quote=F)
