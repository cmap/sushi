library(argparse)
library(tidyverse)

source("utils/kitchen_utensils.R")
source("growth_correction/growth_correction_functions.R")

# Shell script argument parser ----
parser = ArgumentParser()

parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE, help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose", help = "Print little output")
parser$add_argument("--wkdir", default = getwd(), help = "Working directory")
parser$add_argument("--normalized_counts", default = "normalized_counts.csv",
                    help = "Path to file containing normalized counts")
parser$add_argument("--l2fc", default = "l2fc.csv", help = "Path to file containing l2fc values")

args = parser$parse_args()

# init function call from MK
LFC %<>% 
  dplyr::filter(pert_type == "trt_cp") %>% 
  dplyr::group_split(pcr_plate, pcr_well) %>% 
  lapply(apply_growth_correction) %>% 
  dplyr::bind_rows() %>% 
  dplyr::ungroup()