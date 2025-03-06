options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
source("./src/normalize.R")
source("./src/kitchen_utensils.R")
source("./src/qc_table_functions.R")

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
parser$add_argument("--annotated_counts", default="annotated_counts.csv", help = "annotated counts file")
parser$add_argument("--unknown_barcode_counts", default="unknown_barcode_counts.csv", help = "unknown barcode counts file")
parser$add_argument("--cell_set_meta", default="cell_set_and_pool_meta.csv", help = "cell set metadata")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up inputs ----
filtered_counts= data.table::fread(args$filtered_counts, header= TRUE, sep= ',')
CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')
input_pseudocount= as.numeric(args$pseudocount)
input_id_cols= unlist(strsplit(args$id_cols, ","))
annotated_counts = fread(args$annotated_counts, header= TRUE, sep= ',')
unknown_counts = fread(args$unknown_barcode_counts, header= TRUE, sep= ',')
cell_set_meta = fread(args$cell_set_meta, header= TRUE, sep= ',')

# Generate first set of QC flag results
flagged_results = compute_read_stats(filtered_counts= filtered_counts, annotated_counts= annotated_counts,
                                     unknown_counts= unknown_counts, cell_set_meta= cell_set_meta)

# Get dataframe to run through normalization ----
norm_df = flagged_results$result

# Run normalize ----
print("Creating normalized count file ...")
normalized_counts = normalize(X= norm_df, id_cols= input_id_cols,
                              CB_meta= CB_meta, pseudocount= input_pseudocount)

# Write out file ----
normcounts_outpath = paste(args$out, "normalized_counts.csv", sep='/')
print(paste0("Writing normalized count file to ", normcounts_outpath))
write.csv(normalized_counts, normcounts_outpath, row.names= FALSE, quote= FALSE)

# Ensure that normalized file was sucessfully generated ----
check_file_exists(normcounts_outpath)
