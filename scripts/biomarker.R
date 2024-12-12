options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("./src/biomarker_functions.R")
source("./src/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("--collapsed_l2fc", default="l2fc.csv", help="file containing replicate collapsed log fold change values")
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id",
                    help = "Columns that can describe a cell line")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--collapsed_l2fc_column", default="median_l2fc",
                    help = "column containing replicate collapsed log fold change values")
parser$add_argument("--build_dir", default= "", help = "Path to the build directory")
parser$add_argument("--univariate_biomarker", default=true, help="Whether to calculate univariate biomarkers")
parser$add_argument("--multivariate_biomarker", default=true, help="Whether to calculate multivariate biomarkers")
parser$add_argument("--biomarker_file", default="", help="File containing depmap data")

# Get command line options, if help option encountered p3rint help and exit
args <- parser$parse_args()

# Read in l2fc file ----
l2fc= data.table::fread(args$l2fc, header= TRUE, sep= ',')

# Parse some parameters ----
sig_cols= unlist(strsplit(args$sig_cols, ","))
response_column= args$collapsed_l2fc_column
build_dir = args$build_dir
univariate_biomarker = args$univariate_biomarker
multivariate_biomarker = args$multivariate_biomarker
bio_file = args$biomarker_file

# Create treatment_columns by filtering out elements containing "dose"
treatment_columns <- sig_cols[!grepl("dose", sig_cols)]

# Construct output path
out_path <- paste0(build_dir, "/biomarker")

# Run univariate biomarker analysis if requested
if (univariate_biomarker) {
  # Check if the output directory exists, if not create it
  if (!dir.exists(out_path)) {
      dir.create(out_path)
  }
  univariate_table <- create_univariate_biomarker_table(in_path = build_dir,
                                                        out_path = out_path,
                                                        output_file_name = "l2fc_univariate_biomarkers.csv",
                                                        treatment_columns = treatment_columns,
                                                        response_column = response_column,
  )
}

# Run multivariate biomarker analysis if requested
if (multivariate_biomarker) {
  # Check if the output directory exists, if not create it
  if (!dir.exists(out_path)) {
      dir.create(out_path)
  }
  multivariate_table <- create_multivariate_biomarker_table(in_path = build_dir,
                                                            out_path = out_path,
                                                            output_file_name = "l2fc_multivariate_biomarkers.csv",
                                                            treatment_columns = treatment_columns,
                                                            response_column = response_column,
  )
}