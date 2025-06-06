options(cli.unicode = FALSE)
library(argparse)
library(magrittr)
library(tidyverse)
source("biomarker/biomarker_functions.R")
source("utils/kitchen_utensils.R")

# Argument parser ----
parser <- ArgumentParser()
# specify desired options
parser$add_argument("--collapsed_lfc", default="collapsed_l2fc.csv", help="file containing replicate collapsed log fold change values")
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id,lua",
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

# Print some arguments ----
print(paste0("lfc_biomarker is ", lfc_biomarker))
print(paste0("auc_biomarker is " , auc_biomarker))
print(paste0("univariate_biomarker is " , univariate_biomarker))
print(paste0("multivariate_biomarker is " , multivariate_biomarker))

# Create the treatment columns ----
trt_cols = unique(c("day", sig_cols[grepl("pert|is_combination", sig_cols)]))
# ABOVE: Is this reason to include cell set in cell line cols instead of sig_cols?

# Check if the output directory exists, if not create it
if (!dir.exists(out_path)) {
    dir.create(out_path)
}

# Call the biomarker functions ----
if (univariate_biomarker) {
  if (lfc_biomarker) {
    create_univariate_biomarker_table(
      in_path = lfc_path,
      out_path = out_path,
      output_file_name = "median_l2fc_univariate_biomarkers.csv",
      treatment_columns = trt_cols,
      response_column = lfc_column,
      depmap_file = bio_file
    )
  }
  if (auc_biomarker) {
    create_univariate_biomarker_table(
      in_path = drc_file,
      out_path = out_path,
      output_file_name = "log2_auc_univariate_biomarkers.csv",
      treatment_columns = trt_cols,
      response_column = dr_column,
      depmap_file = bio_file
    )
  }
}

if (multivariate_biomarker) {
  if (lfc_biomarker) {
    create_multivariate_biomarker_table(
      in_path = lfc_path,
      out_path = out_path,
      output_file_name = "median_l2fc_multivariate_biomarkers.csv",
      treatment_columns = trt_cols,
      response_column = lfc_column,
      depmap_file = bio_file
    )
  }
  if (auc_biomarker) {
    create_multivariate_biomarker_table(
      in_path = drc_file,
      out_path = out_path,
      output_file_name = "log2_auc_multivariate_biomarkers.csv",
      treatment_columns = trt_cols,
      response_column = dr_column,
      depmap_file = bio_file
    )
  }
}
