#' filter raw reads
#' 
#' takes the raw readcount table and filters for expected indices and cell lines
#' using the given metadata. QC metrics are written out to a text file
#'
#'
#' @param raw_counts - an unfiltered counts table
#' @param sample_meta - the sample metadata for the particular experiment. Must follow the given set of guidelines for metadata
#' @param cell_line_meta - master metadata of cell lines
#' @param cell_set_meta - master metdata of cell sets and their contents
#' @param CB_meta - master metdata of control barcodes, their sequences, and their doses
#' @param id_cols - comma
#' @return list with the following elements
#' #' \itemize{
#'   \item filtered_table df of annotated readcounts
#'   \item qc_table: QC table of index_purity and cell_line_purity 
#' }
#' @export 
filter_raw_reads = function(
  raw_counts, sample_meta, cell_line_meta, 
  cell_set_meta, CB_meta, id_cols=c('treatment', 'dose','dose_unit','day')) {
  
  #insert profile_id into sample_meta
  sample_meta$profile_id = do.call(paste,c(sample_meta[id_cols], sep=':'))
  
  index_filtered = raw_counts %>%
    dplyr::filter(index_1 %in% sample_meta$IndexBarcode1,
                  index_2 %in% sample_meta$IndexBarcode2)
  index_purity = sum(index_filtered$n) / sum(raw_counts$n)

  cell_line_filtered = index_filtered %>%
    merge(sample_meta, by.x=c("index_1", "index_2"), by.y=c("IndexBarcode1", "IndexBarcode2")) %>%
    merge(cell_line_meta, by.x="forward_read_cl_barcode", by.y="Sequence", all.x=T) %>%
    merge(cell_set_meta, by="cell_set", all.x=T) %>%
    dplyr::filter(mapply(grepl, LUA, members) |
                    (mapply(grepl, LUA, cell_set) & is.na(members)) |
                    (forward_read_cl_barcode %in% CB_meta$Sequence))
  cell_line_purity = sum(cell_line_filtered$n) / sum(index_filtered$n)

  qc_table = data.frame(cell_line_purity=cell_line_purity, index_purity = index_purity)
  
  annotated_counts = cell_line_filtered %>%
    merge(CB_meta, by.x="forward_read_cl_barcode", by.y="Sequence", all.x=T) %>%
    dplyr::select_if(function(col) sum(is.na(col)) < length(col)) %>%
    dplyr::select(-any_of(c("flowcell_name", "flowcell_lane", "index_1", "index_2", "members",
                            "lysate_well", "lysate_plate","forward_read_cl_barcode"))) %>%
    dplyr::relocate(project_code, CCLE_name, DepMap_ID, prism_cell_set, Name, log_dose, profile_id, trt_type, control_barcodes,
                    bio_rep, tech_rep) %>%
    dplyr::relocate(n, .after=last_col())
  
  return(list(filtered_counts=annotated_counts, qc_table=qc_table))
}
