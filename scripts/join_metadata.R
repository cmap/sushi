library(argparse)
library(tidyverse)
source("./src/join_metadata.R")

# Argument parser ----
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument('--sample_meta', default= 'sample_meta.csv', help= 'Sample meta data for the sequencing run.')
parser$add_argument("--assay_pool_meta", default="assay_pool_meta.txt", help = "Assay pool metadata")
parser$add_argument('--lfc', default= 'l2fc.csv', help= 'L2FC data.') # level 4
parser$add_argument('--collapsed_lfc', default= 'collapsed_l2fc.csv', help= 'Collapsed l2fc data.') # level 5
parser$add_argument('--sig_cols', default= 'cell_set,treatment,dose,dose_unit,day', 
                    help= 'Columns that uniquely identify a condition.') 
parser$add_argument('--out', default= getwd(), help= 'Path to the output directory.')

args <- parser$parse_args()

# set output to working directory if none is specified
if (args$out == "") {
  args$out = args$wkdir
}

# Read in files and prepare some parameters ----
sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',')
sig_cols= unlist(strsplit(args$sig_cols, ","))

# For assay pool meta, check if it exists. If so, then filter it for relevant cell_sets/davepool_ids
# and select/rename some columns.
assay_pool_meta_exists= FALSE
if(file.exists(args$assay_pool_meta)) {
  assay_pool_meta_exists= TRUE # Update boolean
  
  # Read in assay pool meta and transform the table into something more usable.
  assay_pool_meta= read.delim(args$assay_pool_meta)
  unique_cell_sets= unique(sample_meta$cell_set[sample_meta$cell_set != ""])
  input_assay_pool_meta= assay_pool_meta %>% dplyr::filter(davepool_id %in% unique_cell_sets) %>% 
    dplyr::select(DepMap_ID= depmap_id, CCLE_name= ccle_name, cell_set= davepool_id, pool_id)
}

# Add sample meta and assay pool meta to l2fc table ----
if(file.exists(args$lfc)) {
  print('Attempting to add sample_meta to l2fc file.')
  l2fc= data.table::fread(args$lfc, header= T, sep= ',')
  
  # Add sample meta columns to l2fc
  if('bio_rep' %in% colnames(sample_meta) & 'bio_rep' %in% colnames(l2fc)) {
    input_cols= c(sig_cols, 'bio_rep')
  } else {
    input_cols= sig_cols
    print('WARNING: No "bio_rep" column detected. Proceeding with just sig_cols.')
  }
  l2fc_with_meta_columns= join_metadata(input_df= l2fc, metadata= sample_meta, key_cols= input_cols)
  
  # Add assay pool meta columns to l2fc
  print('Attempting to add assay_pool_meta to l2fc file.')
  if(assay_pool_meta_exists) {
    l2fc_with_meta_columns= join_metadata(input_df= l2fc_with_meta_columns, 
                                          metadata= input_assay_pool_meta,
                                          key_cols= c('DepMap_ID', 'CCLE_name', 'cell_set'))
  } else {
    print('WARNING: Assay pool meta not detected and will not be joined onto l2fc.')
  }
  
  # Write out
  outpath= paste(args$out, 'l2fc_with_meta_columns.csv', sep='/')
  print(paste("Writing l2fc_with_meta_columns.csv to ", outpath))
  l2fc_with_meta_columns %>% write.csv(outpath, row.names= FALSE, quote= FALSE)
} else {
  print('WARNING: l2fc.csv does not exist. Skipping this file.')
}

# Add sample meta and assay pool meta to collapsed_l2fc table ----
if(file.exists(args$collapsed_lfc)) {
  print('Attempting to add sample_meta to collapsed l2fc.')
  collapsed_l2fc= data.table::fread(args$collapsed_lfc, header= T, sep= ',')
  
  # Add sample meta columns to collapsed l2fc
  collapsed_l2fc_with_meta_columns= join_metadata(input_df= collapsed_l2fc, metadata= sample_meta, 
                                                  key_cols= sig_cols)
  
  # Add assay pool meta columns to collapsed l2fc
  if(assay_pool_meta_exists) {
    print('Attempting to add assay_pool_meta to collapsed l2fc.')
    collapsed_l2fc_with_meta_columns= join_metadata(input_df= collapsed_l2fc_with_meta_columns, 
                                                    metadata= input_assay_pool_meta,
                                                    key_cols= c('DepMap_ID', 'CCLE_name', 'cell_set'))
  } else {
    print('WARNING: Assay pool meta not detected and will not be joined onto collapsed l2fc.')
  }
  
  # Write out
  outpath= paste(args$out, 'collapsed_l2fc_with_meta_columns.csv', sep='/')
  print(paste("Writing collapsed_l2fc_with_meta_columns.csv to ", outpath))
  collapsed_l2fc_with_meta_columns %>% write.csv(outpath, row.names= FALSE, quote= FALSE)
} else {
  print('WARNING: collapsed_l2fc.csv does not exist. Skipping this file.')
}
