suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
#suppressPackageStartupMessages(library(scam))
#suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(reshape2))

#TODO fix dependecies


## compute_l2fc
## takes collapsed counts and computes log-fold change values as compared to the annotated control columns
##
## takes:
##      collapsed_counts - table with median_normalized_n column and control_sample column that designates the 
##          name of the control sample for each treatment sample
compute_l2fc = function(collapsed_counts) {
  treatments = collapsed_counts %>% 
    filter(!(trt_type %in% c("trt_ctrl", "negcon", "day_0")))
  controls = collapsed_counts %>% 
    filter(trt_type %in% c("trt_ctrl", "negcon")) %>% 
    mutate(control_sample=sample_ID) %>% 
    rename(control_median_normalized_n = median_normalized_n) %>% 
    dplyr::select(CCLE_name, DepMap_ID, prism_cell_set, control_sample, control_median_normalized_n)
  l2fc = treatments %>% 
    merge(controls, by=c("CCLE_name", "DepMap_ID", "prism_cell_set", "control_sample"), all.x=T, all.y=T) %>% 
    mutate(l2fc=log2(median_normalized_n/control_median_normalized_n)) %>% 
    dplyr::select(CCLE_name, DepMap_ID, prism_cell_set, sample_ID, trt_type, control_sample, 
                  median_normalized_n, control_median_normalized_n, l2fc)
  
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
parser$add_argument("-c", "--collapsed_counts", default="collapsed_counts.csv",
                    help="Path to directory containing fastq files")
parser$add_argument("--out", default="", help = "Output path. Default is working directory")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

collapsed_counts = read.csv(args$collapsed_counts)

print("computing log-fold change")
l2fc = compute_l2fc(collapsed_counts)

write.csv(l2fc, "l2fc.csv", row.names=F, quote=F)