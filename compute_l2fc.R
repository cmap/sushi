suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
#suppressPackageStartupMessages(library(scam))
#suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(reshape2))

## compute_l2fc
## takes normalized counts and computes log-fold change values as compared to the annotated control columns
##
## takes:
##      normalized_counts - table with normalized_n column and control_sample column that designates the 
##          name of the control sample for each treatment sample
compute_l2fc = function(normalized_counts, control_types) {
  treatments = normalized_counts %>% 
    filter(!(trt_type %in% control_types),
           is.na(Name)) %>% 
    dplyr::select(-Name, -log_dose, -n, -log_n, -log_normalized_n) %>% 
    group_by_at(setdiff(names(.), c("normalized_n", "tech_rep", "profile_id"))) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    ungroup()
  controls = normalized_counts %>% 
    filter(trt_type %in% control_types,
           is.na(Name)) %>% 
    mutate(control_sample=sample_ID) %>% 
    dplyr::select(-Name, -log_dose, -n, -log_n, -log_normalized_n) %>% 
    group_by_at(setdiff(names(.), c("normalized_n", "tech_rep", "profile_id"))) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    ungroup() %>% 
    #group_by_at(setdiff(names(.), c("sum_normalized_n", "bio_rep"))) %>% 
    group_by(CCLE_name, DepMap_ID, prism_cell_set, control_sample) %>% 
    dplyr::summarise(control_median_normalized_n = median(sum_normalized_n),
                     control_mad_sqrtN = mad(log10(sum_normalized_n))/sqrt(n())) %>% 
    ungroup() %>% 
    mutate(control_pass_QC = ifelse(control_mad_sqrtN > 0.5, F, T)) %>% 
    dplyr::select(CCLE_name, DepMap_ID, prism_cell_set, control_sample, control_median_normalized_n, control_mad_sqrtN, control_pass_QC)
  l2fc = treatments %>% 
    merge(controls, by=c("CCLE_name", "DepMap_ID", "prism_cell_set", "control_sample"), all.x=T, all.y=T) %>% 
    mutate(l2fc=log2(sum_normalized_n/control_median_normalized_n)) %>% 
    dplyr::relocate(project_code, CCLE_name, DepMap_ID, prism_cell_set, sample_ID, trt_type, control_sample, control_barcodes,
                    bio_rep)
  
  return(l2fc)
}

#Need Arguments
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("-ct", "--control_types", default="trt_ctrl,negcon", help="trt_types to use as control")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

control_types = unlist(strsplit(args$control_types, ","))

normalized_counts = read.csv(args$normalized_counts)

print("computing log-fold change")
l2fc = compute_l2fc(normalized_counts, control_types)

l2fc_out = paste(args$out, "l2fc.csv", sep="/")
write.csv(l2fc, l2fc_out, row.names=F, quote=F)
