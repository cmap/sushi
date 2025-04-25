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
parser$add_argument("--cell_line_cols", default="pool_id,depmap_id",
                    help = "Columns that can describe a cell line")
parser$add_argument("--sig_cols", default="cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("--dose_col", default="pert_dose",
                    help = "column containing dose values")
parser$add_argument("--l2fc_column", default="l2fc",
                    help = "column containing log fold change values")
parser$add_argument("--type_col", default="pert_type",
                    help = "column containing perturbation type, only those with trt_cp will be considered.")
parser$add_argument("--cap_for_viability", default= 1.5, help = "The upper threshold for the viability values for before fitting the curves. Default is 1.5, any viability value above this value will be made equal to 1.5")
parser$add_argument("--build_dir", default= "", help = "Path to the build directory")
parser$add_argument("--out_dir", default= "", help = "Path to the output directory")

# Get command line options, if help option encountered p3rint help and exit
args <- parser$parse_args()

# Read in l2fc file ----
l2fc= data.table::fread(args$l2fc, header= TRUE, sep= ',')

# Parse some parameters ----
cell_line_cols= unlist(strsplit(args$cell_line_cols, ","))
sig_cols= unlist(strsplit(args$sig_cols, ","))
dose_col= args$dose_col
l2fc_column= args$l2fc_column
type_col= args$type_col
cap_for_viability= as.numeric(args$cap_for_viability)
build_dir = args$build_dir
out_dir = args$out_dir

# Create treatment_columns by filtering out elements containing "dose"
treatment_cols <- union(
  sig_cols[grepl("pert", sig_cols) & !grepl("dose", sig_cols)],
  sig_cols[sig_cols == "day"]
)

# Calculate dose response ----
print("Calculating dose response ...")
dose_response= create_drc_table(LFC= l2fc,
                                cell_line_cols= cell_line_cols,
                                treatment_cols= treatment_cols,
                                dose_col= dose_col,
                                l2fc_column= l2fc_column,
                                type_col= type_col,
                                cap_for_viability= cap_for_viability)

# Validation: Check that dose_response is not empty ----
if(nrow(dose_response) == 0) {
  stop('Dose response table is empty.')
}

# Check if the output durectory exists, if not create it
if (!dir.exists(paste0(out_dir))) {
  dir.create(paste0(out_dir))
}

# Write out the DRC table ----
drc_outpath= paste0(out_dir, "/DRC_TABLE.csv")
paste0("Writing DRC_TABLE.csv to ", drc_outpath)
write.csv(dose_response, drc_outpath, row.names= FALSE)

# Check to make sure that the file was generated
check_file_exists(drc_outpath)