#'  QC_images
#'
#'  takes a filtered dataframe of raw read counts and generates and writes out QC images to assess
#'  overall quality of project data
#'
#' @param filtered_counts - dataframe of annotated readcounts that must include the following columns:
#'           n: raw readcounts
#'           profile_id: string unique to each sample as defined by filter_counts method
#'           Name: name of the control barcode that the read corresponds to, or NA (if read is cell line)
#'           CCLE_name: name of the cell line that the read corresponds to, or NA (if read is control barcode)
#'           cell_set: string identifier of cell set expected in a given sample, must match a cell set found in cell_set_meta
#' @param cell_set_meta - a metadata dataframe that contains a mapping from cell set names (e.g. CS5) to 
#'           lists of LUAs in that cell set separated by semicolons
#' @param out - the filepath to the folder in which QC images are meant to be saved, NA by default and 
#'           images are saved in the working directory 
#' @return - NA, QC images are written out to the specified folder
#' @export
QC_images = function(filtered_counts, cell_set_meta, out = NA) {
  if(is.na(out)) {
    out = getwd()
  }
  num_profiles = filtered_counts$profile_id %>% unique() %>% length()
  
  # counts in each sample, colored by cell line or control barcode
  print("generating total_counts image")
  total_counts = filtered_counts %>% ungroup() %>% 
    mutate(type = ifelse(!is.na(CCLE_name), "cell line", "control barcode")) %>% 
    group_by(profile_id, type) %>% 
    dplyr::summarise(total_counts = sum(n))
  
  tc = total_counts %>% 
    ggplot() +
    geom_col(aes(x=profile_id, y=total_counts, fill=type)) +
    labs(x="", y="total counts", fill="") +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=5))
  
  pdf(file=paste(out, "total_counts.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(tc)
  dev.off()
  
  # fraction of expected cell lines in each sample
  print("generating cell_lines_present image")
  num_cell_lines = filtered_counts %>% ungroup() %>% 
    filter(!is.na(CCLE_name)) %>% 
    merge(cell_set_meta, by="cell_set", all.x=T) %>%
    mutate(expected_num_cl = as.character(members) %>% purrr::map(strsplit, ";") %>% purrr::map(`[[`, 1) %>% 
             purrr::map(length) %>% as.numeric()) %>% 
    distinct(CCLE_name, profile_id, cell_set, expected_num_cl) %>% 
    group_by(profile_id, cell_set, expected_num_cl) %>% 
    dplyr::summarise(num_cl = n(),
                     frac_cl = num_cl/expected_num_cl) %>% 
    distinct(profile_id, cell_set, expected_num_cl, num_cl, frac_cl)
  
  ncl = num_cell_lines %>% 
    ggplot() +
    geom_col(aes(x=profile_id, y=frac_cl*100)) +
    labs(x="", y="% expected cell lines present") +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=5))
  
  pdf(file=paste(out, "cell_lines_present.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(ncl)
  dev.off()
  
  # sample correlation
  print("generating sample_cor image")
  correlation_matrix = filtered_counts %>% ungroup() %>% 
    filter(!is.na(CCLE_name)) %>% 
    mutate(log_n = log10(n)) %>% 
    dcast(CCLE_name~profile_id, value.var="log_n") %>% 
    column_to_rownames("CCLE_name") %>% 
    cor(use="pairwise.complete.obs") 
  
  cp = correlation_matrix %>% 
    melt() %>% 
    ggplot() +
    geom_tile(aes(x=Var1, y=Var2, fill=value)) +
    labs(x="", y="", fill="correlation") +
    scale_fill_gradient(low="yellow", high="red") +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=5),
          axis.text.y = element_text(size=5))
  
  pdf(file=paste(out, "sample_cor.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
  print(cp)
  dev.off()
  
  # control barcode trend
  print("generating control_barcode_trend image")
  wells_with_cb = filtered_counts %>% ungroup() %>% 
    filter(control_barcodes=="Y")
  
  if(nrow(wells_with_cb)!=0) {
    cbt = wells_with_cb %>% 
      filter(is.na(CCLE_name)) %>% 
      ggplot(aes(x=log_dose, y=log10(n))) +
      geom_point() +
      geom_smooth(method = 'glm') +
      ggpubr::stat_regline_equation(aes(label =  ..eq.label..), label.y.npc = "top") +
      ggpubr::stat_regline_equation(aes(label =  ..adj.rr.label..), label.y.npc = "top", position = position_nudge(y=-0.5)) +
      facet_wrap(~profile_id) +
      labs(x="log(dose)")
    
    pdf(file=paste(out, "control_barcode_trend.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
    print(cbt)
    dev.off()
  }
  
  # cell line counts vs. control barcode counts
  # assumes that last piece of profile_id is tech_rep
  print("generating all_cell_lines_trend image")
  num_tech_rep = filtered_counts$tech_rep %>% unique() %>% length()
  
  if(num_tech_rep==2) {
    samples_with_two_tech_rep = filtered_counts %>% 
      filter(tech_rep==2) %>% 
      mutate(sample_id = gsub('.{2}$', '', profile_id)) %>% 
      pull(sample_id) 
    samples_with_more_tech_rep = filtered_counts %>% 
      filter(tech_rep>2) %>% 
      mutate(sample_id = gsub('.{2}$', '', profile_id)) %>% 
      pull(sample_id) %>% unique()
    
    if(length(samples_with_more_tech_rep)!=0 & sum(samples_with_more_tech_rep %in% samples_with_two_tech_rep)!=0) {
      samples_with_two_tech_rep = samples_with_two_tech_rep[-which(samples_with_two_tech_rep %in% samples_with_more_tech_rep)]
    }
    
    aclt = filtered_counts %>% ungroup() %>% 
      mutate(type = ifelse(!is.na(CCLE_name), "cell line", "control barcode"),
             log_n = log10(n),
             sample_id = gsub('.{2}$', '', profile_id)) %>% 
      filter(sample_id %in% samples_with_two_tech_rep) %>% 
      dplyr::select(-profile_id, -n) %>% 
      dcast(...~tech_rep, value.var="log_n") %>% 
      ggplot(aes(x=`1`, y=`2`, color=type)) +
      geom_point(alpha=0.5) +
      facet_wrap(~sample_id) +
      labs(x="log10(counts) tech rep 1", y="log10(counts) tech rep 2", color="")
    
    pdf(file=paste(out, "all_cell_lines_trend.pdf", sep="/"),
        width=sqrt(num_profiles), height=sqrt(num_profiles))
    print(aclt)
    dev.off()
  }
}

