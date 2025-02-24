options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("./src/biomarker_functions.R")
source("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("--collapsed_lfc", default="collapsed_l2fc.csv", help="file containing replicate collapsed log fold change values")
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id",
                    help = "Columns that can describe a cell line")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--collapsed_l2fc_column", default="median_l2fc",
                    help = "column containing replicate collapsed log fold change values")
parser$add_argument("--dr_column", default="log2_auc",
                    help = "column containing auc values from dose response")
parser$add_argument("--build_dir", default= "", help = "Path to the build directory")
parser$add_argument("--univariate_biomarker", default="true", help="Whether to calculate univariate biomarkers")
parser$add_argument("--multivariate_biomarker", default="true", help="Whether to calculate multivariate biomarkers")
parser$add_argument("--lfc_biomarker", default="true", help="Whether to calculate lfc biomarkers")
parser$add_argument("--auc_biomarker", default="true", help="Whether to calculate auc biomarkers")
parser$add_argument("--biomarker_file", default="/data/biomarker/current/depmap_datasets_public.h5", help="File containing depmap data")
parser$add_argument("--drc_file", default="DRC_TABLE.csv", help="File containing auc values from dose response")
parser$add_argument("--out_path", default= "", help = "Path to the output directory (per comppound)")

#TODO: Update column headers

# Get command line options, if help option encountered p3rint help and exit
args <- parser$parse_args()

# Basic parameters ----
sig_cols= unlist(strsplit(args$sig_cols, ","))
build_dir = args$build_dir
bio_file = args$biomarker_file

# Paths and columns to use
lfc_path = args$collapsed_lfc
lfc_column= args$collapsed_l2fc_column
dr_column = args$dr_column
drc_file = args$drc_file
out_path = args$out_path

# Parameters for determining which biomarker(s) to calculate
lfc_biomarker = as.logical(toupper(args$lfc_biomarker))
auc_biomarker = as.logical(toupper(args$auc_biomarker))
univariate_biomarker = as.logical(toupper(args$univariate_biomarker))
multivariate_biomarker = as.logical(toupper(args$multivariate_biomarker))

# Create treatment_columns by filtering out elements containing "dose"
treatment_columns <- sig_cols[!grepl("cell_set", sig_cols)]

# Check if the output directory exists, if not create it
if (!dir.exists(out_path)) {
    dir.create(out_path)
}

# Function to determine treatment columns
get_treatment_columns <- function(response_column, sig_cols) {
  if (response_column == args$collapsed_l2fc_column) {
    # Exclude "cell_set" for lfc biomarker
    return(sig_cols[!grepl("cell_set", sig_cols)])
  } else if (response_column == args$dr_column) {
    # Exclude "cell_set" and "dose" for auc biomarker
    return(sig_cols[!grepl("cell_set|dose", sig_cols)])
  } else {
    stop("Unknown response column specified.")
  }
}

# Function to call the creation of biomarker tables
create_biomarker_table <- function(in_path, out_path, response_column, sig_cols, depmap_file, biomarker_type) {
  # Dynamically determine treatment columns
  treatment_columns <- get_treatment_columns(response_column, sig_cols)

  # Choose the appropriate function based on biomarker type
  biomarker_function <- ifelse(biomarker_type == "univariate",
                               create_univariate_biomarker_table,
                               create_multivariate_biomarker_table)

  # Construct the output file name
  output_file_name <- paste0(response_column, "_", biomarker_type, "_biomarkers.csv")

  print(paste0("Creating ", biomarker_type, " biomarker table using ", response_column, " from ", in_path, "..."))

  # Call the appropriate biomarker function
  biomarker_function(
    in_path = in_path,
    out_path = out_path,
    output_file_name = output_file_name,
    treatment_columns = treatment_columns,
    response_column = response_column,
    depmap_file = depmap_file
  )
}

# Process biomarkers based on user inputs
if (univariate_biomarker || multivariate_biomarker) {
  # Pick the datasets and their corresponding response columns
  datasets <- list()
  if (lfc_biomarker) datasets <- c(datasets, list(list(path = lfc_path, response = lfc_column)))
  if (auc_biomarker) datasets <- c(datasets, list(list(path = drc_file, response = dr_column)))

  # Loop through the selected datasets
  for (dataset in datasets) {
    in_path <- dataset$path
    response_column <- dataset$response

    # Run univariate analysis if requested
    if (univariate_biomarker) {
      create_biomarker_table(
        in_path = in_path,
        out_path = out_path,
        response_column = response_column,
        sig_cols = sig_cols,
        depmap_file = bio_file,
        biomarker_type = "univariate"
      )
    }

    # Run multivariate analysis if requested
    if (multivariate_biomarker) {
      create_biomarker_table(
        in_path = in_path,
        out_path = out_path,
        response_column = response_column,
        sig_cols = sig_cols,
        depmap_file = bio_file,
        biomarker_type = "multivariate"
      )
    }
  }
}

