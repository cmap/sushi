options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
source("./src/normalize.R")
source("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--id_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day,bio_rep,tech_rep",
                    help = "Columns to identify each PCR well")
parser$add_argument("--CB_meta", default="CB_meta.csv", help= "Control Barcode metadata")
parser$add_argument("-o", "--out", default=getwd(), help= "Output path. Defaults to working directory")
parser$add_argument("--pseudocount", default=20, help = "pseudocount for normalization")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up inputs ----
filtered_counts= data.table::fread(args$filtered_counts, header= TRUE, sep= ',')
CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')
input_pseudocount= as.numeric(args$pseudocount)
input_id_cols= unlist(strsplit(args$id_cols, ","))

# Run normalize ----
print("Creating normalized count file ...")
normalized_counts = normalize(X= filtered_counts, id_cols= input_id_cols,
                              CB_meta= CB_meta, pseudocount= input_pseudocount)

# Write out file ----
normcounts_outpath = paste(args$out, "normalized_counts.csv", sep='/')
print(paste0("Writing normalized count file to ", normcounts_outpath))
write.csv(normalized_counts, normcounts_outpath, row.names= FALSE, quote= FALSE)

# Ensure that normalized file was sucessfully generated ----
check_file_exists(normcounts_outpath)
