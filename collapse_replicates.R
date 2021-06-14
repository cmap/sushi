suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))


## check_replicate_cor
## checks that technical and biological replicates are all well correlated with each other
## returns normalized counts table, filtered for good replicates, and writes out file reporting bad replicates
## 
## takes:
##      normmalized_counts - dataframe of normalized counts that must included trt_type column 
##          with at least one negcon sample and normalized_n column
check_replicate_cor = function(normalized_counts) {
  tech_rep_cor = normalized_counts %>% 
    filter(trt_type=="negcon", is.na(Name)) %>% 
    dcast(CCLE_name+sample_ID+bio_rep~tech_rep, value.var="normalized_n") %>% 
    dplyr::select(-CCLE_name, -sample_ID, -bio_rep) %>% 
    cor(use="complete.obs") %>% 
    as.data.frame()
  
  write.csv(tech_rep_cor, "tech_rep_cor.csv", row.names=T, quote=F)
  
  tech_collapsed_counts = normalized_counts %>% 
    filter(is.na(Name)) %>%  
    #group_by(CCLE_name, DepMap_ID, prism_cell_set, sample_ID, bio_rep, trt_type, control_sample) %>% 
    group_by(CCLE_name, DepMap_ID, prism_cell_set, sample_ID, bio_rep, trt_type) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    ungroup()
  
  bio_rep_cor = tech_collapsed_counts %>% 
    filter(trt_type=="negcon") %>% 
    dcast(CCLE_name+sample_ID~bio_rep, value.var="sum_normalized_n") %>% 
    dplyr::select(-CCLE_name, -sample_ID) %>% 
    cor(use="complete.obs") %>% 
    as.data.frame()
  
  write.csv(bio_rep_cor, "bio_rep_cor.csv", row.names=T, quote=F)
}

## collapse_counts
## collapses filtered normalized counts and computes MAD/sqrt(n) metrics.
## cell lines with MAD/sqrt(n) > 0.5 are filtered out, and file with filtered out cell lines is written. 
## log10(median counts) vs. MAD/sqrt(n) graph is saved, and collapsed filtered count table is returned
##
## takes:
##      filtered_normalized_counts - normalized counts with bad replicates already filtered out
collapse_counts = function(normalized_counts) {
  collapsed_counts = normalized_counts %>% 
    filter(is.na(Name)) %>%  
    #group_by(CCLE_name, DepMap_ID, prism_cell_set, cell_set, sample_ID, trt_type, control_sample, bio_rep) %>% 
    group_by(CCLE_name, DepMap_ID, prism_cell_set, cell_set, sample_ID, trt_type, bio_rep) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    ungroup() %>% 
    #group_by(CCLE_name, DepMap_ID, prism_cell_set, scell_set, ample_ID, trt_type, control_sample) %>% 
    group_by(CCLE_name, DepMap_ID, prism_cell_set, cell_set, sample_ID, trt_type) %>% 
    dplyr::summarise(median_normalized_n = median(sum_normalized_n),
                     mad_sqrtN = mad(log10(sum_normalized_n))/sqrt(n())) %>% 
    ungroup() %>% 
    mutate(pass_QC = ifelse(mad_sqrtN > 0.5, F, T))
  
  return(collapsed_counts)
}


parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--normalized_counts", default="normalized_counts.csv",
                    help="Path to directory containing fastq files")
parser$add_argument("--out", default="", help = "Output path. Default is working directory")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

normalized_counts = read.csv(args$normalized_counts)

print("checking replicate correlation")
check_replicate_cor(normalized_counts)

print("collapsing counts")
collapsed_counts = collapse_counts(normalized_counts)

collapsed_count_out_file = paste(
  args$out,
  "collapsed_counts.csv",
  sep='/'
)

write.csv(collapsed_counts, collapsed_count_out_file, row.names=F, quote=F)
