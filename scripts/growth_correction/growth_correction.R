library(argparse)
library(tidyverse)

source("utils/kitchen_utensils.R")
source("growth_correction/growth_correction_functions.R")

# Shell script argument parser ----
parser = ArgumentParser()

parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE, help = "Print extra output [default]")
parser$add_argument("-q", "--quietly", action = "store_false", dest = "verbose", help = "Print little output")
parser$add_argument("--l2fc", default = "l2fc.csv", help = "Path to file containing l2fc values.")
parser$add_argument("--growth_annotations", default = "growth_annotations.csv",
                    help = "Path to metadata containing growth annotations.")
parser$add_argument("--sig_cols", default = "cell_set,pert_name,pert_dose,pert_dose_unit,day",
                    help = "columns used to generate signature ids")
parser$add_argument("-o", "--out", default = getwd(), help = "Output path. Default is working directory")

args = parser$parse_args()

# Read in files
growth_annotations = data.table::fread(args$growth_annotations)
l2fc = data.table::fread(args$l2fc)

# Set up some parameters
sig_cols = unlist(strsplit(args$sig_cols, ","))

# Detect if bio_rep column exists - TO DO - MAKE THIS INTO A PARAMETER!
if("bio_rep" %in% colnames(l2fc)) {
  bio_rep_id_cols = c(sig_cols, "bio_rep")
} else {
  bio_rep_id_cols = sig_cols
  message("bio_rep column not detected. Assuming that there are NO biological replicates.")
}

# Add annotations to l2fcs
l2fc[growth_annotations, on = c("pool_id", "lua", "depmap_id"), growth_pattern := i.growth_pattern]

# Throw an error if a cell line does not have a growth annotation.
if (NA %in% unique(l2fc$growth_pattern)) {
  no_annot_cls = unique(l2fc[is.na(growth_pattern), c("pool_id", "lua", "depmap_id")])
  message("The following cell lines are missing growth annotations...")
  print(no_annot_cls)
  stop("Some cell lines do not have growth annotations!")
}

# Correct l2fcs by regressing out cell line growth patterns
corrected_l2fc = l2fc |>
  dplyr::group_split(dplyr::across(tidyrselect::all_of(bio_rep_id_cols))) |>
  lapply(apply_growth_correction, raw_l2fc_col = "l2fc") |>
  dplyr::bind_rows() |>
  dplyr::ungroup() |>
  dplyr::select(-growth_pattern)

# Write out
outpath = file.path(args$out, args$l2fc)
message("Writing corrected l2fc value to ", outpath)
write.csv(corrected_l2fc, outpath, row.names = FALSE, quote = FALSE)