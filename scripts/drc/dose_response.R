options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("drc/dose_response_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("--l2fc", default="l2fc.csv", help="file containing log fold change values")
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id,lua",
                    help = "Columns that can describe a cell line")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--l2fc_column", default="l2fc",
                    help = "column containing log fold change values")
parser$add_argument("--type_col", default="pert_type",
                    help = "column containing perturbation type, only those with trt_cp will be considered.")
parser$add_argument("--cap_for_viability", default = 1.5,
                    help = paste0("The upper threshold for the viability values for before fitting the curves.
                                   Default is 1.5, any viability value above this value will be made equal to 1.5"))
parser$add_argument("--build_dir", default= "", help = "Path to the build directory")
parser$add_argument("--out_dir", default= "", help = "Path to the output directory")

# Get command line options, if help option encountered p3rint help and exit
args <- parser$parse_args()

# Read in l2fc file ----
l2fc= data.table::fread(args$l2fc, header= TRUE, sep= ',')

# Parse some parameters ----
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))
sig_cols= unlist(strsplit(args$sig_cols, ","))
l2fc_col= args$l2fc_col
type_col= args$type_col
cap_for_viability= as.numeric(args$cap_for_viability)
build_dir = args$build_dir

# Pull out dose columns from sig_cols
dose_cols = sig_cols[grepl("_dose$", sig_cols)]

# Print a message if two dose columns were detected
print(paste0("Detecting the following dose column(s): ", dose_cols))
if (length(dose_cols) > 1) {
  print("More than one does column was supplied.")
  # This suggests a CPS run
} else if (length(dose_cols) > 2) {
  stop("More than two dose columns were detected!")
}

# create_drc_table is called for each dose column detected.
# For MTS, the function should just be called once.
# For CPS, the function should be called twice.
drc_outputs = list()

for (idx in seq_along(dose_cols)) {
  print(paste0("Working on ", dose_cols[idx]))

  # Filter l2fc dt by dropping NA in the current dose column
  subset_l2fc = l2fc[!is.na(get(dose_cols[idx])), ]

  # Run drc only if the subsetted dt is not empty. Otherwise continue.
  # This was added for single agents to pass through module post-split
  if (nrow(subset_l2fc) > 0) {
    drc_outputs[[idx]] = create_drc_table(
      # LFC input is filtered only for rows where the current dose column is filled in
      LFC = subset_l2fc,
      cell_line_cols = cell_line_cols,
      treatment_cols = sig_cols[!grepl(dose_cols[idx], sig_cols)],
      dose_col = dose_cols[idx],
      l2fc_col = l2fc_col,
      cap_for_viability = cap_for_viability
    )
  }
}

dose_response = dplyr::bind_rows(drc_outputs)

# Validation: Check that dose_response is not empty ----
if (nrow(dose_response) == 0) {
  warning("Dose response table is empty. An output file will not be created.")
} else {
  # Write out the DRC table ----
  # Check if the output directory exists, if not create it
  if (!dir.exists(args$out_dir)) {
    dir.create(args$out_dir)
  }

  drc_outpath = file.path(args$out_dir, "DRC_TABLE.csv")
  paste0("Writing DRC_TABLE.csv to ", drc_outpath)
  write.csv(dose_response, drc_outpath, row.names = FALSE)

  # Check to make sure that the file was generated
  check_file_exists(drc_outpath)
}
