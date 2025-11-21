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
filtered_counts= data.table::fread(args$filtered_counts, header= TRUE, sep= ',')
CB_meta= data.table::fread(args$CB_meta, header= TRUE, sep= ',')
input_pseudocount= as.numeric(args$pseudocount)
input_id_cols= unlist(strsplit(args$id_cols, ","))
read_detection_limit = as.integer(args$read_detection_limit)
negcon_cols = unlist(strsplit(args$negcon_cols, split = ","))
negcon_type = args$negcon_type
output_file = args$output_file

# Some hard coded inputs
req_negcon_reps = 6

# If normalized counts files already exists, remove them ----
delete_existing_files(args$out, "normalized_counts")

# Hotfix: Catch dev runs without negative controls
if (pseudocount == 0 & !negcon_type %in% unique(filtered_counts$pert_type)) {
  message("Error: The pseudocount is set to zero.")
  message("If this is a NextSeq run change the pseudocount from 0 to 20.")
  stop("Pseudocount parameter in the config needs to be non zero if there are no negative controls!")
}

# Run normalize ----
print("Creating normalized count file ...")
normalized_counts = normalize(X = filtered_counts, id_cols = input_id_cols, req_negcon_reps = req_negcon_reps,
                              negcon_type = negcon_type, CB_meta = CB_meta, pseudocount = input_pseudocount)

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
    # Error out if no pseudocount is provided, but there are not enough negative control replicates
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
write.csv(normalized_counts, normcounts_outpath, row.names= FALSE, quote= FALSE)

# Ensure that normalized file was sucessfully generated ----
check_file_exists(normcounts_outpath)
