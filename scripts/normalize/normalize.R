options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
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
parser$add_argument("--read_detection_limit", default = "10",
                    help = "Smallest read count value. Used to compute pseudovalue.")
parser$add_argument("--negcon_cols", default = "pcr_plate,pert_vehicle",
                    help = "List of columns in filtered counts that describe a negative control condition.")
parser$add_argument("--negcon_type", default = "ctl_vehicle",
                    help = "String in the column pert_type that indicates a negative control.")
parser$add_argument("--pseudocount", default = 0, help = "Pseudocount used in normalization.")
parser$add_argument("--output_file", default = "normalized_counts.csv", help = "File name for normalized counts.")
parser$add_argument("-o", "--out", default=getwd(), help= "Output path. Defaults to working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# Set up inputs ----
filtered_counts= read_data_table(args$filtered_counts)
CB_meta= read_data_table(args$CB_meta)
input_pseudocount= as.numeric(args$pseudocount)
input_id_cols= unlist(strsplit(args$id_cols, ","))
read_detection_limit = as.integer(args$read_detection_limit)
negcon_cols = unlist(strsplit(args$negcon_cols, split = ","))
negcon_type = args$negcon_type
output_file = args$output_file

# Hard coded input -
# This is the number of negative control replicates required to perform
# the MAD filter on the control barcodes and to add a plate specific pseudocount.
req_negcon_reps = 6

# Remove normalized files that already exist in the directory
delete_existing_files(args$out, "normalized_counts")

# Identify CBs in PCR wells that can be normalized.
message("Identifying control barcodes to normalize ...")
cbs_to_norm = get_valid_norm_cbs(filtered_count = filtered_counts,
                                 CB_meta = CB_meta,
                                 id_cols = input_id_cols,
                                 negcon_type = negcon_type,
                                 cb_mad_cutoff = 1,
                                 req_negcon_reps = req_negcon_reps)

normalized_counts = normalize(X = filtered_counts,
                              valid_cbs = cbs_to_norm,
                              id_cols = input_id_cols,
                              pseudocount = input_pseudocount)

# Check if pseudovalue addition is needed
if (input_pseudocount < read_detection_limit) {
  # Determine the number of negative control replicates
  negcon_reps = filtered_counts[pert_type == negcon_type] |>
    dplyr::distinct(across(all_of(input_id_cols))) |>
    dplyr::group_by(pcr_plate) |>
    dplyr::summarise(num_negcon_reps = dplyr::n(), .groups = "drop")

  if (min(negcon_reps$num_negcon_reps) >= req_negcon_reps) {
    # Run pseudovalue addition if enough negcon reps are present
    message("Adding pseudovalue corresponding to a read count of ", read_detection_limit, "...")
    normalized_counts = add_pseudovalue(normalized_counts, negcon_cols = negcon_cols,
                                        read_detection_limit = read_detection_limit,
                                        negcon_type = negcon_type)
  } else {
    # Error out if no pseudocount is provided and there are not enough negative control replicates
    message("Not enough negative control replicates for pseudovalue calculations.")
    message("Provide a pseudocount for normalization instead. ")
    stop("Not enough negative control replicates for every PCR plate.")
  }

} else {
  # Continue without pseudovalue addition if a large pseudocount is provided.
  message("Read counts normalized with a pseudocount of ", input_pseudocount, ".")
  message("Pseudovalue addition was skipped due to high pseudocount provided for normalization.")
}

# Write out file ----
normcounts_outpath = file.path(args$out, output_file)
print(paste0("Writing normalized count file to ", normcounts_outpath))
write_out_table(normalized_counts, normcounts_outpath)

# Ensure that normalized file was sucessfully generated ----
check_file_exists(normcounts_outpath)