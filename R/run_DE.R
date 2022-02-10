#'  run_DE
#'
#'  takes a filtered dataframe of raw read counts and applies DEseq to generate
#'  logfold change, standard errors, and significance values
#'
#' @param filtered_counts - dataframe of annotated readcounts that must include the following columns:
#'           n: raw readcounts
#'           trt_type: spcification of sample type (e.g. whether a sample is a control or not)
#'           CCLE_name: contains the name of the cell line that the read corresponds to, or NA
#'           Name: contains the name of the control barcode that the read corresponds to, or NA
#'           all columns specified in sample_cols and sig_cols 
#' @param sample_cols - a vector of column names denoting which values specify each individual sample
#'                    (by default: cell_set, treatment, dose, dose_unity, day, and bio_rep)
#' @param sig_cols - a vector of column names denoting which values specify each individual signature
#'                    (by default: cell_set, treatment, dose, dose_unity, and day)
#' @param control_type - character value specifying which samples to use as controls. Must match a value in the 
#'                    trt_type column
#' @return Table with logfold change values, standard errors, and significance values for all signatures
#' @export
run_DE = function(filtered_counts, 
                  sample_cols = c("cell_set", "treatment", "dose", "dose_unit", "day", "bio_rep"), 
                  sig_cols = c("cell_set", "treatment", "dose", "dose_unit", "day"), 
                  control_type = "negcon") {
  filtered_counts$sample_id = do.call(paste,c(filtered_counts[sample_cols], sep=':'))
  filtered_counts$sig_id = do.call(paste,c(filtered_counts[sig_cols], sep=':'))
  
  countData = filtered_counts %>% 
    filter(!trt_type %in% c("day_0", "empty")) %>% 
    dplyr::select(-any_of(c("log_dose", "log_n"))) %>% 
    dplyr::group_by_at(setdiff(names(.), c("n", "tech_rep", "profile_id", "pcr_plate", "pcr_well"))) %>% 
    dplyr::summarise(n=sum(n)) %>% ungroup() %>% 
    mutate(CCLE_name = ifelse(is.na(CCLE_name), as.character(Name), as.character(CCLE_name)),
           CCLE_name = paste0(CCLE_name, "__", cell_set)) %>% 
    dcast(CCLE_name~sample_id, value.var="n") %>% 
    column_to_rownames("CCLE_name") %>% 
    as.matrix()
  countData[is.na(countData)] = 0
  
  colData = filtered_counts %>% 
    filter(!trt_type %in% c("day_0", "empty")) %>% 
    mutate(sig_id = ifelse(trt_type==control_type, "control", as.character(sig_id))) %>% 
    distinct(sample_id, sig_id) %>% 
    arrange(match(sample_id, colnames(countData))) %>% 
    column_to_rownames("sample_id") %>% 
    mutate(sig_id = relevel(as.factor(sig_id), ref="control"))
  
  BCs = rownames(countData) %>% grepl("^BC", .) %>% unlist()
  
  dds = DESeqDataSetFromMatrix(countData = countData,
                               colData = colData,
                               design = ~sig_id)
  dds = estimateSizeFactors(dds, controlGenes = BCs)
  dds = DESeq(dds)
  
  final = data.frame()
  for(sig_id in unique(colData$sig_id)){
    if(sig_id!="control") {
      comparison = paste0("sig_id_", str_replace_all(sig_id, regex("[[:punct:]-[_]+[\\s+]]"), "."), "_vs_control") #"-|:|%"
      print(comparison)
      #res_shrunk = lfcShrink(dds, coef=comparison, type="apeglm")
      res_shrunk = results(dds, name=comparison)
      hold = res_shrunk[,c("log2FoldChange", "lfcSE", "padj")] %>% 
        as.data.frame() %>% 
        mutate(sig_id = sig_id) %>% 
        rownames_to_column("CCLE_name")
      final = final %>% 
        rbind(hold)
    }
  }
  final = final %>% 
    mutate(cell_set = as.character(CCLE_name) %>% purrr::map(str_split, "__") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 2) %>% unlist() %>% as.character(),
           CCLE_name = as.character(CCLE_name) %>% purrr::map(str_split, "__") %>% purrr::map(`[[`, 1) %>% purrr::map(`[`, 1) %>% unlist() %>% as.character())
  
  l2fc_DEseq = filtered_counts %>% 
    filter(!trt_type %in% c("day_0", "empty", control_type)) %>% 
    mutate(CCLE_name = ifelse(is.na(CCLE_name), as.character(Name), as.character(CCLE_name))) %>%
    select(-any_of(c("Name", "log_dose", "n", "log_n", "bio_rep", "tech_rep", "profile_id", "sample_id"))) %>% 
    distinct() %>% 
    merge(final, by=c("CCLE_name", "cell_set", "sig_id"), all.x=F)
  
  return(l2fc_DEseq)
}

  
  