suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))

#source("../src/load_libraries.R")
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr)) #write_delim
suppressPackageStartupMessages(library(stringr)) #str_detect
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(tidyr)) #pivot_wider
suppressPackageStartupMessages(library(reshape2))
library(prismSeqR)


## print_args
## writes configuration to file
##
## takes:
##      args: args object from argparse
print_args <- function(args){
  config <- data.frame(args=names(args), values=unname(unlist(args)))
  config_path = paste(
    args$out,
    "config.txt",
    sep="/"
  )
  print(paste("Saving config.txt file in :", config_path))
  write_delim(config, config_path, delim = ": ", col_names=F)
}

# create parser object
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                    dest="verbose", help="Print little output")

parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--raw_counts", default="raw_counts.csv", help = "path to file containing raw counts")
parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")
parser$add_argument("-s", "--sample_meta", default="", help = "Sample metadata")
parser$add_argument("--cell_line_meta", default="../metadata/cell_line_meta.csv", help = "Cell Line metadata")
parser$add_argument("--cell_set_meta", default="../metadata/cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")
parser$add_argument("--id_cols", default="cell_set,treatment,dose,dose_unit,day,bio_rep,tech_rep",
    help = "Columns used to generate profile ids, comma-separated colnames from --sample_meta")
parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
parser$add_argument("--reverse_index2", action="store_true", default=FALSE, help = "Reverse complement of index 2 for NovaSeq")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

#Save arguments
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

CB_meta = read.csv(args$CB_meta)
cell_line_meta = read.csv(args$cell_line_meta)
cell_set_meta = read.csv(args$cell_set_meta)
sample_meta = read.csv(args$sample_meta)
raw_counts = read.csv(args$raw_counts)

#split id_cols args
id_cols = unlist(strsplit(args$id_cols, ","))

if (!all(id_cols %in% colnames(sample_meta))){
  stop(paste("Colnames not found in sample_meta, check metadata or --id_cols argument:", args$id_cols))
}

sample_meta$profile_id = do.call(paste,c(sample_meta[id_cols], sep=':'))

print("creating filtered count file")
filtered_counts = filter_raw_reads(
  raw_counts,
  sample_meta,
  cell_line_meta,
  cell_set_meta,
  CB_meta,
  id_cols=id_cols,
  count_threshold= count_threshold,
  reverse_index2=args$reverse_index2
)

# Write out module outputs
qc_table = filtered_counts$qc_table
qc_out_file = paste(args$out, 'QC_table.csv', sep='/')
print(paste("writing QC_table to: ", qc_out_file))
write.csv(qc_table, qc_out_file, row.names=F, quote=F)

filtered_counts = filtered_counts$filtered_counts
filtrc_out_file = paste(args$out, 'filtered_counts.csv', sep='/')
print(paste("writing filtered counts csv to: ", filtrc_out_file))
write.csv(filtered_counts, filtrc_out_file, row.names=F, quote=F)

annotated_counts = filtered_counts$annotated_counts
annot_out_file = paste(args$out, 'annotated_counts.csv', sep='/')
print(paste("writing annotated counts to: ", annot_out_file))
write.csv(annotated_counts, annot_out_file, row.names=F, quote=F)
