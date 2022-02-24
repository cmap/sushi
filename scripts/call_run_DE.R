suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(prismSeqR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))

parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--sample_cols", default="cell_set,treatment,dose,dose_unit,day,bio_rep", 
                    help = "columns used to generate sample ids")
parser$add_argument("--sig_cols", default="cell_set,treatment,dose,dose_unit,day", 
                    help = "columns used to generate signature ids")
parser$add_argument("-ct", "--control_type", default="negcon", help="trt_type to use as control")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

filtered_counts = read.csv(args$filtered_counts)
sample_cols = unlist(strsplit(args$sample_cols, ","))
sig_cols = unlist(strsplit(args$sig_cols, ","))
control_type = args$control_type

if (!all(sample_cols %in% colnames(filtered_counts))){
  stop(paste("Colnames not found in filtered_counts, check metadata or --sample_cols argument:", 
             paste(args$sample_cols, collapse=", "), "\n"))
}
if (!all(sig_cols %in% colnames(filtered_counts))){
  stop(paste("Colnames not found in filtered_counts, check metadata or --sig_cols argument:", 
             paste(args$sig_cols, collapse=", "), "\n"))
}
if (!(control_type %in% filtered_counts$trt_type)){
  stop(paste("control type not found in filtered_counts, check metadata or --control_type argument:", 
             args$control_type, "\n"))
}

print("running DESeq on filtered counts")

DE_out = data.frame()
for(cs in unique(filtered_counts$cell_set)) {
  subset = filtered_counts %>% 
    filter(cell_set==cs)
  hold = run_DE(subset,
                sample_cols,
                sig_cols,
                control_type)
  DE_out = DE_out %>% 
    rbind(hold)
}

DE_out_file = paste(
  args$out,
  "DESeq_l2fc.csv",
  sep='/'
)

if (args$verbose){
  print(paste("Writing DESeq l2fc file", DE_out_file))
}

write.csv(DE_out, DE_out_file, row.names=F, quote=F)


