options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
source("normalize/normalize_functions.R")
source("utils/kitchen_utensils.R")
source("qc_tables/qc_tables_functions.R")

# Argument parser ----
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--id_cols", default="pcr_plate,pcr_well",
                    help = "Columns to identify each PCR well")
parser$add_argument("--CB_meta", default="CB_meta.csv", help= "Control Barcode metadata")
parser$add_argument("--min_read_count", default = "10", 
                    help = "Smallest read count value. Used to compute pseudovalue.")
parser$add_argument("--negcon_cols", default = "pcr_plate,pert_vehicle",
                    help = "List of columns in filtered counts that describe a negative control condition.")
parser$add_argument("--negcon_type", default = "ctl_vehicle",
                    help = "String in the column pert_type that indicates a negative control.")
parser$add_argument("--pseudocount", default = 20, help = "Pseudocount used in normalization.")
parser$add_argument("--output_file", default = "normalized_counts.csv", help = "File name for normalized counts.")
parser$add_argument("-o", "--out", default=getwd(), help= "Output path. Defaults to working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up inputs ----
filtered_counts= data.table::fread(args$filtered_counts, header= TRUE, sep= ',')
CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')
input_pseudocount= as.numeric(args$pseudocount)
input_id_cols= unlist(strsplit(args$id_cols, ","))
min_read_count = as.integer(args$min_read_count)
negcon_cols = unlist(strsplit(args$negcon_cols, split = ","))
negcon_type = args$negcon_type
output_file = args$output_file

# Run normalize ----
print("Creating normalized count file ...")
normalized_counts = normalize(X= filtered_counts, id_cols= input_id_cols,
                              CB_meta= CB_meta, pseudocount= input_pseudocount)
# Add pseudovalue
message("Adding pseudovalue corresponding to a read count of ", min_read_count, "...")
normalized_counts = add_pseudovalue(normalized_counts, negcon_cols = negcon_cols,
                                    min_read_count = min_read_count, negcon_type = negcon_type)

# Write out file ----
normcounts_outpath = file.path(args$out, output_file)
print(paste0("Writing normalized count file to ", normcounts_outpath))
write.csv(normalized_counts, normcounts_outpath, row.names= FALSE, quote= FALSE)

# Ensure that normalized file was sucessfully generated ----
check_file_exists(normcounts_outpath)
