suppressPackageStartupMessages(library(argparse))
suppressMessages(library(cmapR))
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
    cor(use="complete.obs")
  
  bad_tech_rep = tech_rep_cor %>% 
    as.data.frame() %>% 
    rownames_to_column("tech_rep") %>% 
    melt(id.vars="tech_rep") %>% 
    filter(tech_rep!=variable) %>% 
    group_by(tech_rep) %>% 
    dplyr::summarise(median_cor = median(value)) %>% 
    filter(median_cor < 0.8) %>% 
    pull(tech_rep) %>% 
    as.numeric()
  
  if(length(bad_tech_rep)!=0) {
    filtered_normalized_counts = normalized_counts %>% 
      filter(!(tech_rep %in% bad_tech_rep))
  } else {
    filtered_normalized_counts = normalized_counts
  }
  
  bio_rep_cor = filtered_normalized_counts %>% 
    filter(trt_type=="negcon", is.na(Name)) %>% 
    dcast(CCLE_name+sample_ID+tech_rep~bio_rep, value.var="normalized_n") %>% 
    dplyr::select(-CCLE_name, -sample_ID, -tech_rep) %>% 
    cor(use="complete.obs")
  
  bad_bio_rep = bio_rep_cor %>% 
    as.data.frame() %>% 
    rownames_to_column("bio_rep") %>% 
    melt(id.vars="bio_rep") %>% 
    filter(bio_rep!=variable) %>% 
    group_by(bio_rep) %>% 
    dplyr::summarise(median_cor = median(value)) %>% 
    filter(median_cor < 0.8) %>% 
    pull(bio_rep) %>% 
    as.numeric()
  
  if(length(bad_bio_rep)!=0) {
    filtered_normalized_counts = filtered_normalized_counts %>% 
      filter(!(bio_rep %in% bad_bio_rep))
  } else {
    filtered_normalized_counts = filtered_normalized_counts
  }
  
  if(length(bad_tech_rep)==0 & length(bad_bio_rep)==0) {
    write("tech reps removed: None \nbio reps removed: None",
          "replicate_QC.txt")
  } else {
    write(paste0("tech reps removed: ", bad_tech_rep, "\n",
                 "bio_reps_removed: ", bad_bio_rep),
          "replicate_QC.txt")
  }
  
  return(filtered_normalized_counts)
}

## collapse_counts
## collapses filtered normalized counts and computes MAD/sqrt(n) metrics.
## cell lines with MAD/sqrt(n) > 0.5 are filtered out, and file with filtered out cell lines is written. 
## log10(median counts) vs. MAD/sqrt(n) graph is saved, and collapsed filtered count table is returned
##
## takes:
##      filtered_normalized_counts - normalized counts with bad replicates already filtered out
collapse_counts = function(filtered_normalized_counts) {
  collapsed_counts = filtered_normalized_counts %>% 
    filter(is.na(Names)) %>%  
    group_by(CCLE_name, DepMap_ID, prism_cell_set, sample_ID, trt_type, control_sample) %>% 
    dplyr::summarise(median_normalized_n = median(normalized_n),
                     mad_sqrtN = mad(log10(normalized_n))/sqrt(n())) %>% 
    ungroup()
  
  bad_cell_lines = collapsed_counts %>% 
    filter(trt_type=="negcon") %>% 
    filter(mad_sqrtN > 0.5) %>% 
    pull(CCLE_name) %>% 
    as.character()
  
  if(length(bad_cell_lines)==0) {
    write("cell lines removed: None",
          "cell_line_QC.txt")
  } else {
    write(paste0("cell lines removed: ", bad_cell_lines),
          "cell_line_QC.txt")
  }
  
  p = collapsed_counts %>% 
    filter(trt_type=="negcon") %>% 
    ggplot() +
    geom_point(aes(x=log10(median_normalized_n), y=mad_sqrtN)) +
    labs(x="log10(median normalized counts)", y="MAD/sqrt(n)")
  ggsave("median_v_MADsqrtN.jpeg", p, jpeg())
  
  filtered_collapsed_counts = collapsed_counts %>% 
    filter(!(CCLE_name %in% bad_cell_lines)) %>% 
    dplyr::select(-mad_sqrtN)
  
  return(filtered_collapsed_counts)
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
filtered_normalized_counts = check_replicate_cor(normalized_counts)

print("collapsing counts")
collapsed_counts = collapse_counts(filtered_normalized_counts)

collapsed_count_out_file = paste(
  args$out,
  "collapsed_counts.csv",
  sep='/'
)

write.csv(collapsed_counts, collapsed_count_out_file, row.names=F, quote=F)
