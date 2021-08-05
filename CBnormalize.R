suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
## normalize
## takes a filtered dataframe of raw read counts and normalizes
## counts using control barcodes
##
## takes:
##      X - dataframe of annotated readcounts that must include the following columns:
##          log_n: log10 of read counts
##          log_dose: log10 of dose at which control barcode was spiked in, if applicable
##          sample_ID: some identifier that distinguishes between each sample
##          Name: contains the name of the control barcode that the read corresponds to, or NA
##      barcodes - a vector of control barcode Name identifiers
normalize <- function(X, barcodes) {
  normalized <- X %>%
    dplyr::group_by(sample_ID) %>%
    dplyr::mutate(log_normalized_n = glm(y ~ x,
                                         data = tibble(
                                           y = log_dose[Name %in% barcodes],
                                           x = log_n[Name %in% barcodes])) %>% 
                    predict(newdata = tibble(x = log_n))) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(normalized_n = 10^log_normalized_n)
  
  return(normalized)
}

parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")


parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")



# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

filtered_counts = read.csv(args$filtered_counts)
CB_meta = read.csv(args$CB_meta)

print("creating normalized count file")
normalized_counts = filtered_counts %>% 
  mutate(log_n = log10(n)) %>% 
  normalize(CB_meta$Name) 

normcounts_out_file = paste(
  args$out,
  "normalized_counts.csv",
  sep='/'
)

if (args$verbose){
  print(paste("Writing normalized count file", normcounts_out_file))
}

write.csv(normalized_counts, normcounts_out_file, row.names=F, quote=F)


