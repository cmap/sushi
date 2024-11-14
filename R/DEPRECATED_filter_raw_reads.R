#' filter raw reads
#' 
#' takes the raw readcount table and filters for expected indices and cell lines
#' using the given metadata. QC metrics are returned as a data.frame
#'
#' @param raw_counts - an unfiltered counts table
#' @param sample_meta - the sample metadata for the particular experiment. Must follow the given set of 
#'                      guidelines for metadata
#' @param cell_line_meta - master metadata of cell lines
#' @param cell_set_meta - master metdata of cell sets and their contents
#' @param CB_meta - master metdata of control barcodes, their sequences, and their doses
#' @param id_cols - vector of column names used to generate unique profile_id for each sample. 
#'                  cell_set,treatment,dose,dose_unit,day,bio_rep,tech_rep by default
#' @return - list with the following elements
#' #' \itemize{
#'   \item filtered_table df of annotated readcounts
#'   \item qc_table: QC table of index_purity and cell_line_purity 
#' }
#' @export 
filter_raw_reads = function(
  raw_counts, sample_meta, cell_line_meta, 
  cell_set_meta, CB_meta,
  id_cols=c('cell_set','treatment','dose','dose_unit','day','bio_rep','tech_rep'),
  reverse_index2= FALSE, count_threshold= 40) {
  
  # New: convert CB_meta from log10 to log2
  CB_meta= CB_meta %>% dplyr::mutate(log2_dose= log_dose/log10(2)) %>% dplyr::select(-log_dose)
  print("Converting CB_meta from log10 to log2")
  
  if (reverse_index2) {
    sample_meta$IndexBarcode2 <- chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$IndexBarcode2))
    print("Reverse-complementing index 2 barcode.")
  }
  
  print("Filtering raw counts")
  index_filtered = raw_counts %>%
    dplyr::filter(index_1 %in% sample_meta$IndexBarcode1, index_2 %in% sample_meta$IndexBarcode2)
  print("Computing index purity")
  index_purity = sum(index_filtered$n) / sum(raw_counts$n)
  
  print("Filtering cell lines")
  cell_line_filtered = index_filtered %>%
    merge(sample_meta, by.x=c("index_1", "index_2"), by.y=c("IndexBarcode1", "IndexBarcode2")) %>%
    merge(cell_line_meta, by.x="forward_read_barcode", by.y="dna_sequence", all.x=T) %>% # NEW
    merge(cell_set_meta, by="cell_set", all.x=T) %>%
    dplyr::filter(mapply(grepl, lua, members) | # NEW
                    (mapply(grepl, lua, cell_set) & is.na(members)) | # NEW
                    (forward_read_barcode %in% CB_meta$Sequence))
  cell_line_purity = sum(cell_line_filtered$n)/ sum(index_filtered$n)
  
  print("Generating QC table ...")
  qc_table = data.frame(cell_line_purity=cell_line_purity, index_purity = index_purity)
  
  # make template of expected reads
  #index_to_well= sample_meta %>% dplyr::distinct(pick(c('IndexBarcode1', 'IndexBarcode2', 'pcr_plate', 'pcr_well')))
  sample_meta$profile_id= do.call(paste,c(sample_meta[id_cols], sep=':'))
  
  template= sample_meta %>% merge(cell_set_meta, by='cell_set', all.x=T) %>%
    dplyr::mutate(members= ifelse(is.na(members), str_split(cell_set, ';'), str_split(members, ';'))) %>% 
    unnest(cols=c(members)) %>%
    merge(cell_line_meta, by.x= 'members', by.y= 'lua', all.x= T) # NEW
  
  # check for control barcodes and add them to the template
  if ('Y' %in% sample_meta$control_barcodes | T %in% sample_meta$control_barcodes) {
    cb_template= sample_meta %>% dplyr::filter(control_barcodes %in% c('Y', 'T', T)) %>%
      dplyr::mutate(joiner= 'temp') %>%
      merge(CB_meta %>% dplyr::mutate(joiner= 'temp'), by='joiner') %>% dplyr::select(-joiner)
    template= plyr::rbind.fill(template, cb_template)
  }
  
  # annotating reads now takes much longer
  print("Annotating reads ...")
  annotated_counts= raw_counts %>%
    merge(cell_line_meta, by.x="forward_read_barcode", by.y="dna_sequence", all.x=T) %>% # NEW
    merge(CB_meta, by.x="forward_read_barcode", by.y="Sequence", all.x=T) %>%
    merge(sample_meta, by.x= c('index_1', 'index_2'), by.y= c('IndexBarcode1', 'IndexBarcode2'), all.x=T) %>%
    merge(template %>% dplyr::mutate(expected_read= T), 
          by.x= c('index_1', 'index_2', 'forward_read_barcode', intersect(colnames(template), colnames(.))),
          by.y= c('IndexBarcode1', 'IndexBarcode2', 'dna_sequence', intersect(colnames(template), colnames(.))), # NEW
          all.x=T, all.y=T) %>% 
    dplyr::mutate(n= replace_na(n, 0),
                  expected_read= replace_na(expected_read, F))
  
  # filtered counts
  print("Filtering reads ...")
  filt_cols= c('project_code', 'pcr_plate', 'pcr_well', 'ccle_name', 'depmap_id', 'prism_cell_set',
               'control_barcodes', 'Name', 'log2_dose','profile_id', 'trt_type')
  filtered_counts= annotated_counts %>% dplyr::filter(expected_read) %>%
    dplyr::select(any_of(c(filt_cols, id_cols, 'n'))) %>%
    dplyr::mutate(flag= ifelse(n==0, 'Missing', NA),
                  flag= ifelse(n!=0 & n < count_threshold, 'low counts', flag))
  
  # Adjusted column naming for case sensitivity - may need to adjust across several modules
  filtered_counts <- filtered_counts %>%
    rename("CCLE_name" = "ccle_name",
           "DepMap_ID" = "depmap_id")

  annotated_counts <- annotated_counts %>%
    rename("CCLE_name" = "ccle_name",
           "DepMap_ID" = "depmap_id")
  # excluded counts
  #excluded_counts= annotated_counts %>% dplyr::filter(is.na(project_code)) %>%
  #  dplyr::select_if(function(col) sum(is.na(col)) < length(col)) # ignore columns with all NAs
  
  return(list(annotated_counts= annotated_counts, filtered_counts= filtered_counts,
              qc_table= qc_table))
}
