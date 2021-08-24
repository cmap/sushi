suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(cdsrbiomarker))

## collapse_counts
## collapses filtered normalized counts and computes MAD/sqrt(n) metrics.
## cell lines with MAD/sqrt(n) > 0.5 are filtered out, and file with filtered out cell lines is written. 
## log10(median counts) vs. MAD/sqrt(n) graph is saved, and collapsed filtered count table is returned
##
## takes:
##      filtered_normalized_counts - normalized counts with bad replicates already filtered out
generate_biomarkers = function(collapsed_values) {
  bio_in = collapsed_values %>% 
    filter(trt_pass_QC) %>% 
    dcast(DepMap_ID~sprofile_id, value.var="median_l2fc") %>% 
    column_to_rownames("DepMap_ID")
  
  bio_out = cdsrbiomarker::get_biomarkers(bio_in)
  
  return(bio_out)
}


parser <- ArgumentParser()
# specify our desired options 
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
