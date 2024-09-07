options(cli.unicode = FALSE)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr)) #write_delim
suppressPackageStartupMessages(library(stringr)) #str_detect
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(tidyr)) #pivot_wider
suppressPackageStartupMessages(library(sets))
suppressPackageStartupMessages(library(tidyverse)) # load last - after dplyr
source("./src/filter_raw_reads.R")
filtered_counts_required_cols <- c("CCLE_name", "cell_set", "DepMap_ID")
assay_pool_required_cols <- c("ccle_name", "davepool_id", "depmap_id", "pool_id")
cell_line_required_cols <- c("CCLE_name", "DepMap_ID", "LUA", "Sequence")

##### FUNCTIONS #####

#' validate_seq_index_cols
#' 
#' This function checks that a list of sequencing index columns are present in the sample_meta.
#' 
#' @param sequencing_index_cols A vector of strings each representing a column name
#' @param sample_meta Sample meta as a dataframe
validate_seq_index_cols = function(sequencing_index_cols, sample_meta){
  sequencing_index_cols_split = unlist(strsplit(sequencing_index_cols, ","))
  if (!all(sequencing_index_cols_split %in% colnames(sample_meta))){
    stop(paste("All seq columns not found in sample_meta, check metadata or --sequencing_index_cols argument:",
               sequencing_index_cols_split))
  }
}

# TODO
validate_cell_line_meta = function(cell_line_meta) {
  validate_required_cols(cell_line_meta, cell_line_required_cols)
  # stop(paste("Cell line meta validation not yet implemented."))
}

remove_duplicate_luas = function(cell_line_meta) {
  # make sure LUA codes in cell line meta are unique
  cell_line_meta %<>% 
    dplyr::group_by(LUA) %>% 
    dplyr::mutate(LUA.duplicity = n()) %>% 
    dplyr::ungroup()
  
  print(paste0("LUAs that are duplicated ", 
               dplyr::filter(cell_line_meta, LUA.duplicity > 1)$LUA %>% 
                 unique() %>% sort() %>% paste0(collapse = ", "))) # print LUA duplicates
  
  cell_line_meta %<>% 
    dplyr::filter(!duplicated(cell_line_meta$LUA, fromLast = TRUE)) %>%
    dplyr::select(-LUA.duplicity)
  
  return(cell_line_meta)
}

merge_pool_info = function(df, assay_pool_meta, unique_cell_sets){
  assay_pool_subset <- assay_pool_meta[assay_pool_meta$davepool_id %in% unique_cell_sets,] %>% 
    select(pool_id, ccle_name, davepool_id, depmap_id)
  result = df %>% merge(assay_pool_subset, by.x=filtered_counts_required_cols, 
                         by.y=assay_pool_required_cols, all.x=T)
  return(result)
  
}

validate_required_cols = function(df, required_columns){

  if (length(intersect(colnames(df), required_columns)) != length(required_columns)){
    stop(paste("Required columns:", required_columns, "are not present in the provided dataframe."))
  }
}

remove_data_from_filter_counts = function(filtered_counts, outpath) {
  # Remove data if needed
  print('rm_data is TRUE, removing supplied data.')
  data_to_remove <- read.csv(paste(outpath, 'data_to_remove.csv', sep='/'))
  print('Data to remove:')
  print(head(data_to_remove))
  filt_rm <- remove_data(filtered_counts, data_to_remove)
  return(filt_rm)
}

write_results_to_csv = function(results, outpath) {
  # Write out module outputs ----
  qc_table = results$qc_table
  qc_out_file = paste(outpath, 'QC_table.csv', sep='/')
  print(paste("writing QC_table to: ", qc_out_file))
  write.csv(qc_table, qc_out_file, row.names=F, quote=F)
  
  unmapped_reads= results$unmapped_reads
  unmapped_out = paste(outpath, 'unmapped_reads.csv', sep='/')
  print(paste("writing unmapped reads to: ", unmapped_out))
  write.csv(unmapped_reads, unmapped_out, row.names=F)
  
  annotated_counts = results$annotated_counts
  annot_out_file = paste(outpath, 'annotated_counts.csv', sep='/')
  print(paste("writing annotated counts to: ", annot_out_file))
  write.csv(annotated_counts, annot_out_file, row.names=F)
  
  filtered_counts = results$filtered_counts
  filtrc_out_file = paste(outpath, 'filtered_counts.csv', sep='/')
  print(paste("writing filtered counts csv to: ", filtrc_out_file))
  write.csv(filtered_counts, filtrc_out_file, row.names=F, quote=F)
  
  if (results$filtered_counts_original) {
    # keep the full filtered counts with the data that needs to be removed
    write.csv(results$filtered_counts_original, paste(outpath, 'filtered_counts_original.csv', sep='/'), row.names=F, quote=F)
  }
}
readFiles = function(args){
  cell_set_meta= data.table::fread(args$cell_set_meta, header= T, sep= ',', data.table= F)
  cell_line_meta= data.table::fread(args$cell_line_meta, header= T, sep= ',', data.table= F)
  CB_meta= data.table::fread(args$CB_meta, header= T, sep= ',', data.table= F)
  sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',', data.table= F)
  raw_counts= data.table::fread(args$raw_counts, header= T, sep= ',', data.table= F)

  //IS THIS CORRECT
  ds <- data.frame(cell_set_meta, cell_line_meta,CB_meta,sample_meta,raw_counts)

  return ds
  
}
parseArgs = function(){
   parser <- ArgumentParser()
    # specify desired options
  parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                      help="Print extra output [default]")
  parser$add_argument("-q", "--quietly", action="store_false",
                      dest="verbose", help="Print little output")
  parser$add_argument("--wkdir", default=getwd(), help="Working directory")
  parser$add_argument("-c", "--raw_counts", default="raw_counts.csv", help = "path to file containing raw counts")
  parser$add_argument("-o", "--out", default="", help = "Output path. Default is working directory")
  parser$add_argument("-s", "--sample_meta", default="sample_meta.csv", help = "Sample metadata")
  parser$add_argument("--cell_line_meta", default="cell_line_meta.csv", help = "Cell Line metadata")
  parser$add_argument("--cell_set_meta", default="cell_set_meta.csv", help = "Cell set metadata")
  parser$add_argument("--assay_pool_meta", default="assay_pool_meta.txt", help = "Assay pool metadata")
  parser$add_argument("--CB_meta", default="CB_meta.csv", help = "Control Barcode metadata")
  parser$add_argument("--sequencing_index_cols", default= "index_1,index_2", 
                      help = "Sequencing columns in the sample meta")
  parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
  parser$add_argument("--reverse_index2", type="logical",
                      help = "Reverse complement of index 2 for NovaSeq and NextSeq")
  parser$add_argument("--rm_data", type="logical", help = "Remove bad experimental data")
  parser$add_argument("--pool_id", type="logical", help = "Pull pool IDs from CellDB.")
  parser$add_argument("--control_type", default="negcon", 
                      help = "negative control wells in trt_type column in sample metadata")
  
  # get command line options, if help option encountered print help and exit
  return parser$parse_args()
}
main = function() {
  
  # get command line options, if help option encountered print help and exit
  args <- parseArgs()
  
  # set output to working directory if none is specified
  if (args$out == ""){
    args$out = args$wkdir
  }
  #print_args(args)
  
  # Read in files and set up parameters ----
  #cell_set_meta= data.table::fread(args$cell_set_meta, header= T, sep= ',', data.table= F)
  #cell_line_meta= data.table::fread(args$cell_line_meta, header= T, sep= ',', data.table= F)
  #CB_meta= data.table::fread(args$CB_meta, header= T, sep= ',', data.table= F)
  #sample_meta= data.table::fread(args$sample_meta, header= T, sep= ',', data.table= F)
  #raw_counts= data.table::fread(args$raw_counts, header= T, sep= ',', data.table= F)

  dfc = readFiles(args);
  cell_set_meta= dfc$cell_set_meta
  cell_line_meta= dfc$cell_line_meta
  CB_meta= dfc$CB_meta
  sample_meta= dfc$sample_meta
  raw_counts= dfc$raw_counts

  validate_cell_line_meta(cell_line_meta)
  validate_seq_index_cols(args$sequencing_index_cols, sample_meta)
  cell_line_meta = remove_duplicate_luas(cell_line_meta)
  
  # Run filter_raw_reads -----
  print("creating filtered count file")
  filtered_counts = filter_raw_reads(raw_counts,
                                     sample_meta,
                                     cell_line_meta,
                                     cell_set_meta,
                                     CB_meta,
                                     sequencing_index_cols= sequencing_index_cols,
                                     count_threshold= as.numeric(args$count_threshold),
                                     reverse_index2= args$reverse_index2)
  
  validate_required_cols(filtered_counts$filtered_counts, filtered_counts_required_cols)
  unique_cell_sets <- unique(sample_meta$cell_set[sample_meta$cell_set != ""])
  
  if (args$pool_id) {
    assay_pool_meta = read.delim(args$assay_pool_meta)
    validate_required_cols(assay_pool_meta, assay_pool_required_cols)
    filtered_counts$filtered_counts = merge_pool_info(filtered_counts$filtered_counts, assay_pool_meta, unique_cell_sets)
    filtered_counts$annotated_counts = merge_pool_info(filtered_counts$annotated_counts, assay_pool_meta, unique_cell_sets)
  }
  
  
  # Validation: Basic file size check ----
  if(sum(filtered_counts$filtered_counts$n) == 0) {
    stop('All entries in filtered counts are missing!')
  }
  
  cl_entries= filtered_counts$filtered_counts %>% dplyr::filter(!is.na(CCLE_name))
  if(sum(cl_entries$n) == 0) {
    stop('All cell line counts are zero!')
  }
  
  print(paste("rm_data:", args$rm_data))
  if(args$rm_data == TRUE){
    filtered_counts$filtered_counts_original <- filtered_counts$filtered_counts
    filtered_counts$filtered_counts = remove_data_from_filter_counts(filtered_counts$filtered_counts, args$out)
    
    # re-point to what filtered_counts should be
    filtered_counts <- filt_rm
    rows_removed = nrow(filtered_counts_original) - nrow(filtered_counts)
    paste("Number of rows removed: ", rows_removed)
  }
  
  write_results_to_csv(filtered_counts, args$out)
}

# main()
