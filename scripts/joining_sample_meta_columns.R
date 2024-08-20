library(argparse)
library(tidyverse)
source("./src/join_sample_meta.R")

# Argument parser ----
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument('--sample_meta', default= 'sample_meta.csv', help= 'Sample meta data for the sequencing run.')
parser$add_argument('--l2fc', default= 'l2fc.csv', help= 'L2FC data.') # level 4
parser$add_argument('--collapsed_l2fc', default= 'collapsed_l2fc.csv', help= 'Collapsed l2fc data.') # level 5
parser$add_argument('--sig_cols', default= 'cell_set,treatment,dose,dose_unit,day', 
                    help= 'Columns that uniquely identify a condition.') 
parser$add_argument('--out', default= getwd(), help= 'Path to the output directory.')

args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == "") {
  args$out = args$wkdir
}

# Prepare args ----
sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',')
sig_cols= unlist(strsplit(args$sig_cols, ","))

# Add in metadata for l2fc file ----
if(file.exists(args$l2fc)) {
  l2fc= data.table::fread(args$l2fc, header= T, sep= ',')
  if('bio_rep' %in% sample_meta & 'bio_rep' %in% l2fc) {
    input_cols= c(sig_cols, 'bio_rep')
  } else {
    input_cols= sig_cols
    print('WARNING: No "bio_rep" column detected. Proceeding with just sig_cols.')
  }
  l2fc_with_sm= join_sample_meta(df= l2fc, sample_meta, key_cols= input_cols)
  
  # Write out
  outpath= paste(args$out, 'l2fc_with_sm.csv', sep='/')
  print(paste("Writing l2fc_with_sm.csv to ", outpath))
  write.csv(l2fc_with_sm, outpath, row.names= FALSE, quote= FALSE)
} else {
  print('WARNING: l2fc.csv does not exist. Skipping this file.')
}
#

# Add in metadata for collapsed_l2fc file ----
if(file.exists(args$collapsed_l2fc)) {
  collapsed_l2fc= data.table::fread(args$collapsed_l2fc, header= T, sep= ',')
  collapsed_l2fc_with_sm= join_sample_meta(df= collapsed_l2fc, sample_meta, key_cols= sig_cols)
  
  # Write out
  outpath= paste(args$out, 'collapsed_l2fc_with_sm.csv', sep='/')
  print(paste("Writing collapsed_l2fc_with_sm.csv to ", outpath))
  write.csv(collapsed_l2fc_with_sm, outpath, row.names= FALSE, quote= FALSE)
} else {
  print('WARNING: collapsed_l2fc.csv does not exist. Skipping this file.')
}
#
