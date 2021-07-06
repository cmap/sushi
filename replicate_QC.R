suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tibble))

## check_replicate_cor
## checks that technical and biological replicates are all well correlated with each other
## returns normalized counts table, filtered for good replicates, and writes out file reporting bad replicates
## 
## takes:
##      normmalized_counts - dataframe of normalized counts that must included trt_type column 
##          with at least one negcon sample and normalized_n column
check_replicate_cor = function(normalized_counts, out) {
  tech_rep_cor = normalized_counts %>% 
    filter(is.na(Name)) %>% 
    dcast(CCLE_name~sample_ID+bio_rep+tech_rep, value.var="log_normalized_n") %>% 
    dplyr::select(-CCLE_name) %>% 
    cor(use="complete.obs") %>% as.data.frame() 
  
  trep_out = paste(args$out, "tech_rep_cor.csv", sep='/')
  write.csv(tech_rep_cor, trep_out, row.names=T, quote=F)
  
  tech_rep_cor_long = tech_rep_cor %>% 
    rownames_to_column("sample_1") %>% 
    melt(id.vars="sample_1", variable.name="sample_2", value.name="cor") %>% 
    mutate(sample_ID_1 = as.character(sample_1) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 1) %>% unlist(),
           sample_ID_2 = as.character(sample_2) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 1) %>% unlist()) %>% 
    filter(sample_ID_1 == sample_ID_2) %>% 
    mutate(bio_rep_1 = as.character(sample_1) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 2) %>% unlist(),
           bio_rep_2 = as.character(sample_2) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 2) %>% unlist()) %>% 
    filter(bio_rep_1 == bio_rep_2) %>%
    mutate(tech_rep_1 = as.character(sample_1) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 3) %>% unlist(),
           tech_rep_2 = as.character(sample_2) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 3) %>% unlist()) %>%
    filter(tech_rep_2>tech_rep_1) %>% 
    dplyr::rename(sample_ID = sample_ID_1, bio_rep = bio_rep_1) %>% 
    dcast(sample_ID+bio_rep~tech_rep_1+tech_rep_2, value.var="cor")
  
  trep_long_out = paste(args$out, "tech_rep_cor_long.csv", sep='/')
  write.csv(tech_rep_cor_long, trep_long_out, row.names=T, quote=F)
  
  tech_collapsed_counts = normalized_counts %>% 
    filter(is.na(Name)) %>%  
    dplyr::select(-Name, -log_dose, -n, -log_n, -log_normalized_n, -profile_id) %>% 
    group_by_at(setdiff(names(.), c("normalized_n", "tech_rep"))) %>% 
    dplyr::summarise(sum_normalized_n = sum(normalized_n)) %>% 
    ungroup()
  
  bio_rep_cor = tech_collapsed_counts %>% 
    dcast(CCLE_name~sample_ID+bio_rep, value.var="sum_normalized_n") %>% 
    dplyr::select(-CCLE_name) %>% 
    cor(use="complete.obs") %>% 
    as.data.frame()
  
  brep_out = paste(args$out, "bio_rep_cor.csv", sep='/')
  write.csv(bio_rep_cor, brep_out, row.names=T, quote=F)
  
  bio_rep_cor_long = bio_rep_cor %>% 
    rownames_to_column("sample_1") %>% 
    melt(id.vars="sample_1", variable.name="sample_2", value.name="cor") %>% 
    mutate(sample_ID_1 = as.character(sample_1) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 1) %>% unlist(),
           sample_ID_2 = as.character(sample_2) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 1) %>% unlist()) %>% 
    filter(sample_ID_1 == sample_ID_2) %>% 
    mutate(bio_rep_1 = as.character(sample_1) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 2) %>% unlist(),
           bio_rep_2 = as.character(sample_2) %>% purrr::map(strsplit, "_") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 2) %>% unlist()) %>% 
    filter(bio_rep_2>bio_rep_1) %>% 
    dplyr::rename(sample_ID = sample_ID_1) %>% 
    dcast(sample_ID~bio_rep_1+bio_rep_2, value.var="cor")
  
  brep_long_out = paste(args$out, "bio_rep_cor_long.csv", sep='/')
  write.csv(bio_rep_cor_long, brep_long_out, row.names=T, quote=F)
}


parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

normalized_counts = read.csv(args$normalized_counts)

print("checking replicate correlation")
check_replicate_cor(normalized_counts, args$out)
