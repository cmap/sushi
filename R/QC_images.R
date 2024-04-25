#'  QC_images
#'
#'  takes a filtered dataframe of raw read counts and generates and writes out QC images to assess
#'  overall quality of project data
#'  
#' @param sample_meta - sample metadata
#' @param annotated_counts -
#' @param filtered_counts - dataframe of annotated readcounts that must include the following columns:
#'           n: raw readcounts
#'           profile_id: string unique to each sample as defined by filter_counts method
#'           Name: name of the control barcode that the read corresponds to, or NA (if read is cell line)
#'           CCLE_name: name of the cell line that the read corresponds to, or NA (if read is control barcode)
#'           cell_set: string identifier of cell set expected in a given sample, must match a cell set 
#'           found in cell_set_meta
#' @param normalized_counts - 
#' @param CB_meta - control barcode metadata
#' @param cell_set_meta - a metadata dataframe that contains a mapping from cell set names (e.g. CS5) to 
#'           lists of LUAs in that cell set separated by semicolons
#' @param out - the filepath to the folder in which QC images are meant to be saved, NA by default and 
#'           images are saved in the working directory 
#' @param sig_cols - 
#' @param count_col_names -
#' @param count_threshold - threshold for low counts
#' @param reverse_index2 reverse index 2 if newer sequencers are used.
#' @return - NA, QC images are written out to the specified folder
#' @export

QC_images = function(sample_meta, annotated_counts, filtered_counts, normalized_counts= NA,
                     CB_meta, cell_set_meta, out = NA, sig_cols, count_col_name= 'normalized_n',
                     count_threshold= 40,
                     reverse_index2= FALSE) {
  if(is.na(out)) {
    out = getwd()
  }
  num_profiles = filtered_counts$profile_id %>% unique() %>% length()
  
  # Some preprocessing ----
  # Reverse index 2 barcodes
  if(reverse_index2) {
    print("Reverse-complementing index 2 barcode.")
    sample_meta$IndexBarcode2= chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$IndexBarcode2))
  }
  
  cb_check= filtered_counts %>%
    dplyr::filter(control_barcodes %in% c("Y", "T", T),
                  !(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type))
  contains_cbs= ifelse(nrow(cb_check)!= 0, T, F)
  #
  
  # Sequencing QCs ____________________ ----
  ## Index count summaries ----
  print("generating index counts tables")
  # pull unique indices.
  expected_index1= unique(sample_meta$IndexBarcode1)
  expected_index2= unique(sample_meta$IndexBarcode2)
  
  index1_counts= annotated_counts %>% dplyr::group_by(index_1) %>%
    dplyr::summarise(idx_n= sum(n, na.rm= T)) %>% dplyr::ungroup() %>%
    dplyr::mutate(fraction= idx_n/sum(idx_n),
                  expected= ifelse(index_1 %in% expected_index1, T, F),
                  contains_n= ifelse(grepl('N', index_1), T, F),
                  lv_dist= apply(stringdist::stringdistmatrix(index_1, expected_index1, method="lv"), 1, min),
                  ham_dist= apply(stringdist::stringdistmatrix(index_1, expected_index1, method="hamming"), 1, min)) %>%
    dplyr::arrange(desc(fraction))
  index2_counts= annotated_counts %>% dplyr::group_by(index_2) %>%
    dplyr::summarise(n= sum(idx_n, na.rm= T)) %>% dplyr::ungroup() %>%
    dplyr::mutate(fraction= idx_n/sum(idx_n),
                  expected= ifelse(index_2 %in% expected_index2, T, F),
                  contains_n= ifelse(grepl('N', index_2), T, F),
                  lv_dist= apply(stringdist::stringdistmatrix(index_2, expected_index2, method="lv"), 1, min),
                  ham_dist= apply(stringdist::stringdistmatrix(index_2, expected_index2, method="hamming"), 1, min)) %>%
    dplyr::arrange(desc(fraction))
  
  # export
  index1_counts %>% write.csv(file= paste(out, 'index1_counts.csv', sep='/'), row.names=F)
  index2_counts %>% write.csv(file= paste(out, 'index2_counts.csv', sep='/'), row.names=F)
  
  ## Total counts ----
  print("generating total_counts image")
  total_counts= filtered_counts %>% ungroup() %>% 
    mutate(type = ifelse(!is.na(CCLE_name), "cell line", "control barcode"),
           sample_id= paste(pcr_plate, pcr_well, sep='_')) %>%
    group_by(pcr_plate, pcr_well, sample_id, profile_id, type) %>% 
    dplyr::summarise(total_counts = sum(n))
  
  tc= total_counts %>% ggplot() +
    geom_col(aes(x=sample_id, y=total_counts, fill=type), alpha=0.75, position='identity') +
    geom_hline(yintercept= 10^4, linetype=2) + 
    facet_wrap(~pcr_plate, scale= 'free_x') +
    labs(x="", y="total counts", fill="") + theme_bw() +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=5)) 
  
  pdf(file=paste(out, "total_counts.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(tc)
  dev.off()
  
  ## Control barcode counts ----
  if(contains_cbs) { 
    print('generating control barcode counts')
    idx_pairs_cbs= annotated_counts %>% dplyr::filter(!is.na(Name), expected_read) %>%
      dplyr::group_by(index_1, index_2, pcr_plate, pcr_well) %>%
      dplyr::summarise(total_well_cbs= sum(n)) %>%
      dplyr::mutate(rpm= total_well_cbs/10^6,
                    p_row= str_extract(pcr_well, '^[a-z,A-Z]'),
                    p_col= str_extract(pcr_well,'\\d+') %>% as.numeric()) %>% ungroup() %>%
      dplyr::arrange(pcr_plate, p_row, p_col) %>% dplyr::select(-p_row, -p_col)
    
    idx_pairs_cbs %>% write.csv(file= paste(out, 'cbs_counts.csv', sep='/'), row.names=F)
  }
  
  ## Distribution of included and excluded reads ----
  distrib_fc= annotated_counts %>% dplyr::filter(n > 1, !is.na(pcr_plate), !is.na(pcr_well)) %>%
    ggplot(aes(x= log10(n), fill= expected_read)) +
    geom_histogram(position='identity', alpha=0.5) +
    geom_vline(xintercept= log10(count_threshold), linetype=2) +
    facet_wrap(~pcr_plate, scales='free') +
    labs(x= 'log10(n)', title= 'Distribution of excluded and included reads') +
    theme_bw()
  pdf(file=paste(out, "count_distribution.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(distrib_fc)
  dev.off()
  
  ## In set/out set reads by well ----
  print('generating in set/out set figures')
  idx_pairs= annotated_counts %>% dplyr::mutate(in_set= ifelse(expected_read, T, F)) %>%
    dplyr::group_by(index_1, index_2, pcr_plate, pcr_well, in_set) %>% 
    dplyr::summarise(sum_n= sum(n, na.rm=T)) %>% dplyr::ungroup() %>%
    dplyr::group_by(index_1, index_2) %>% 
    dplyr::mutate(fraction= sum_n/ sum(sum_n),
                  rpm= sum_n/10^6,
                  p_row= substr(pcr_well, 1, 1),
                  p_col= str_extract(pcr_well, '\\d+') %>% as.numeric()) %>% ungroup()
  
  in_set_percent= idx_pairs %>% dplyr::filter(!is.na(pcr_well), in_set) %>%
    ggplot(aes(x= as.factor(p_col), y=reorder(p_row, desc(p_row)), fill=fraction, label= round(fraction, 2))) + 
    geom_tile(color= 'white') + geom_text(size=3) +
    scale_fill_gradientn(breaks= c(0.5, 0.75, 1), 
                         colours= c('white','yellow','red'),
                         limits=c(0.5,1), oob=squish)+
    coord_fixed(1) + facet_wrap(pcr_plate~.) +
    labs(x= 'plate column', y= 'plate row', title= 'Fraction of \'in set\' reads') +
    theme_bw() #+ theme(legend.position = "bottom")
  
  pdf(file=paste(out, "in_set_percent.pdf", sep="/"), 
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(in_set_percent)
  dev.off()
  
  out_set_rpm= idx_pairs %>% dplyr::filter(!is.na(pcr_well), !in_set) %>%
    ggplot(aes(x= as.factor(p_col), y=reorder(p_row, desc(p_row)), fill= rpm, label= round(rpm, 2))) + 
    geom_tile(color='white') + geom_text(size=3) +
    scale_fill_gradient2(low= 'white', high='red') +
    coord_fixed(1) + facet_wrap(pcr_plate~.) +
    labs(x= 'plate column', y= 'plate row', title= 'Out of set reads (RPM)') +
    theme_bw() #+ theme(legend.position= 'bottom')
  
  pdf(file=paste(out, "out_set_rpm.pdf", sep="/"), 
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(out_set_rpm)
  dev.off()
  
  # Assay QCs _________________________ ----
  ## Cell lines recovered ----
  print("generating cell_lines_present image")
  recovery= filtered_counts %>% ungroup() %>% 
    dplyr::filter(!is.na(CCLE_name)) %>%
    merge(cell_set_meta, by="cell_set", all.x=T) %>%
    dplyr::mutate(members= ifelse(is.na(members), cell_set, members), # for custom cell sets
                  expected_num_cl = as.character(members) %>% purrr::map(strsplit, ";") %>% 
                    purrr::map(`[[`, 1) %>% purrr::map(length) %>% as.numeric(),
                  count_type= ifelse(n > count_threshold, 'Detected', 'Low'),
                  count_type= ifelse(n==0, 'Missing', count_type)) %>%
    group_by(profile_id, pcr_plate, pcr_well, count_type, expected_num_cl) %>% 
    dplyr::summarise(count= n()) %>% 
    dplyr::mutate(frac_type= count/expected_num_cl)
  
  cl_rec= recovery %>% ggplot() +
    geom_col(aes(x=profile_id, y=frac_type*100, fill= reorder(count_type, desc(count_type)))) +
    facet_wrap(~pcr_plate, scales= 'free_x') +
    labs(x="", y="Percentage of expected cell lines", fill= '') +
    theme_bw() +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=5))
  
  pdf(file=paste(out, "cell_lines_present.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(cl_rec)
  dev.off()
  
  ## Fraction of cell lines with low counts ----
  low_cls= filtered_counts %>% dplyr::filter(n < count_threshold, is.na(Name)) %>%
    dplyr::group_by(project_code, pcr_plate, pcr_well, cell_set) %>%
    dplyr::summarize(num_low_cl= length(unique(DepMap_ID))) %>% dplyr::ungroup() %>%
    merge(cell_set_meta, by="cell_set", all.x=T) %>%
    dplyr::mutate(members= ifelse(is.na(members), cell_set, members), # for custom cell sets
                  expected_num_cl = as.character(members) %>% purrr::map(strsplit, ";") %>% 
                    purrr::map(`[[`, 1) %>% purrr::map(length) %>% as.numeric(),
                  fraction= num_low_cl/expected_num_cl,
                  p_row= substr(pcr_well, 1, 1),
                  p_col= str_extract(pcr_well, '\\d+') %>% as.numeric())
  
  low_cl_percent= low_cls %>% ggplot(aes(x= as.factor(p_col), y= reorder(p_row, desc(p_row)), 
                                         fill= fraction, label= round(fraction,2))) +
    geom_tile(color ='white') + geom_text(size=3) +
    scale_fill_gradientn(breaks= c(0, 1), 
                         colours= c('white', 'red'),
                         limits= c(0, 1), oob=squish) +
    coord_fixed(1) + facet_wrap(pcr_plate~.) +
    labs(x= 'plate column', y= 'plate row', fill= 'Fraction of \nlow lines',
         title= 'Fraction of cell lines with low counts') +
    theme_bw()
  
  pdf(file=paste(out, "low_cl_percent.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(low_cl_percent)
  dev.off()
  
  recurring_low_cl= filtered_counts %>% dplyr::filter(n < count_threshold, is.na(Name), trt_type=='negcon') %>%
    dplyr::group_by(project_code, CCLE_name, DepMap_ID, cell_set) %>%
    dplyr::summarize(num_negcon_profiles= n(),
                     median_n= median(n), max_n= max(n)) %>% dplyr::ungroup() %>%
    dplyr::arrange(desc(num_negcon_profiles))
  recurring_low_cl %>% write.csv(file= paste(out, 'recurring_low_cl.csv', sep='/'), row.names=F)
  #
  ## Contaminants ----
  contams= annotated_counts %>% 
    dplyr::filter(!is.na(pcr_plate), !is.na(pcr_well), expected_read==F, !is.na(CCLE_name) | !is.na(Name)) %>%
    dplyr::mutate(barcode_id= ifelse(is.na(CCLE_name), Name, CCLE_name)) %>%
    dplyr::group_by(forward_read_cl_barcode, barcode_id) %>% 
    dplyr::summarise(num_wells= n(), median_n=median(n), max_n= max(n)) %>% ungroup() %>%
    dplyr::arrange(desc(num_wells))
  contams %>% write.csv(file= paste(out, 'contams.csv', sep='/'), row.names=F)
  #
  
  ## Cumulative counts by lines in negcons ----
  print("generating cummulative image")
  cdf= filtered_counts %>% ungroup() %>% 
    dplyr::filter(trt_type=='negcon') %>% # filter for only negative controls
    merge(cell_set_meta, by="cell_set", all.x=T) %>%
    dplyr::mutate(members= ifelse(is.na(members), cell_set, members), # for custom cell cets
                  expected_num_cl= as.character(members) %>% purrr::map(strsplit, ';') %>% 
                    purrr::map(`[[`, 1) %>% purrr::map(length) %>% as.numeric(),
                  expected_num_cl= ifelse(contains_cbs, expected_num_cl + length(unique(CB_meta$Name)),
                                          expected_num_cl)) %>% # add CBs to expected_num_cl if there are CBs
    dplyr::group_by(pcr_plate, pcr_well, profile_id, expected_num_cl) %>% 
    dplyr::mutate(total_counts= sum(n), pct_counts= n/total_counts,) %>% dplyr::arrange(-n) %>% 
    dplyr::mutate(cum_pct= cumsum(pct_counts), rank= row_number(),
                  rank_pct= rank/expected_num_cl) %>% ungroup() %>%
    dplyr::select(DepMap_ID, Name, log2_dose, pcr_plate, pcr_well, profile_id, expected_num_cl, pct_counts,
                  cum_pct, rank, rank_pct)
  
  # additional tables
  mark50= cdf %>% dplyr::filter(cum_pct >= 0.5) %>% dplyr::group_by(profile_id) %>% 
    arrange(cum_pct) %>% dplyr::filter(row_number()==1) %>% ungroup() %>% 
    dplyr::select(profile_id, rank_pct= rank_pct, num50= rank, num50_loc= rank_pct)
  mark95= cdf %>% dplyr::group_by(profile_id) %>% 
    dplyr::mutate(auc= sum(cum_pct*(1/expected_num_cl))) %>% # calculate AUCs
    dplyr::filter(cum_pct >= 0.95) %>% 
    arrange(cum_pct) %>% dplyr::filter(row_number() ==1) %>% ungroup() %>% 
    dplyr::select(profile_id, rank_pct= rank_pct, num95= rank, num95_loc= rank_pct, auc)
  
  cdf_plot= cdf %>% 
    merge(mark50, by= c('profile_id', 'rank_pct'), all.x=T) %>% 
    merge(mark95, by= c('profile_id', 'rank_pct'), all.x= T) %>%
    ggplot(aes(x= rank_pct, y=cum_pct)) +
    { if(contains_cbs) geom_point(. %>% dplyr::filter(!is.na(Name)), 
               mapping=aes(x= rank_pct, y=cum_pct, color=reorder(Name, log2_dose)), size=3) } + 
    geom_line(color='black') +
    # point for 50% of counts
    geom_segment(aes(x= -Inf , y= .50, xend = num50_loc, yend = .50), color= 'black', linetype='dashed') +
    geom_segment(aes(x= num50_loc, y= -Inf, xend = num50_loc, yend = .50), color= 'black', linetype='dashed') +
    geom_label(aes(x=num50_loc, y= .25, label= num50), hjust= 0, color= 'black') +
    # point for 95% of counts
    geom_segment(aes(x= -Inf , y= .95, xend = num95_loc, yend = .95), color= 'black', linetype='dashed') +
    geom_segment(aes(x= num95_loc, y= -Inf, xend = num95_loc, yend = .95), color= 'black', linetype='dashed') +
    geom_label(aes(x=num95_loc, y= .75, label= num95), hjust= 0, color= 'black') +
    # label for AUC
    geom_label(aes(x=num95_loc, y= .25, label= paste0('AUC ', round(auc,3))), hjust= 'inward', color= 'black') +
    facet_wrap(~profile_id) + 
    labs(x='% rank of unique reads', y='Cumulative percentage', color= 'CBs') + theme_bw()
  
  pdf(file=paste(out, "cdf_plot.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(cdf_plot)
  dev.off()
  
  ## Control barcode trends ----
  if(contains_cbs & is.data.frame(normalized_counts)) {
    print("generating control_barcode_trend image")
    cb_trend= normalized_counts %>% 
      dplyr::filter(control_barcodes %in% c("Y", "T", T),
                    !(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type),
                    !is.na(Name)) %>%
      dplyr::group_by(profile_id) %>%
      dplyr::mutate(mean_y= mean(log2_dose),
                    residual2= (log2_dose - log2_normalized_n)^2,
                    squares2= (log2_dose - mean_y)^2,
                    r2= 1 - sum(residual2)/sum(squares2)) %>% ungroup()
    
    trend_sc= cb_trend %>% dplyr::mutate(profile_id= reorder(profile_id, r2)) %>%
      ggplot(aes(x=log2_n, y=log2_dose)) + geom_point() +
      geom_abline(aes(slope=1, intercept= cb_intercept) , color='blue') +
      geom_text(aes(x= min(log2_n), y= desc(sort(unique(log2_dose)))[1], label= paste('r2=', round(r2, 4), sep='')), 
                hjust='inward', vjust='inward') +
      facet_wrap(~profile_id, scales= 'free_x') +
      labs(x= 'log2(n)', y= 'log2(dose)') + theme_bw()
    
    pdf(file=paste(out, "control_barcode_trend.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
    print(trend_sc)
    dev.off()
    
    trend_hist= cb_trend %>% 
      dplyr::distinct(project_code, profile_id, pcr_plate, pcr_well, cb_intercept, r2) %>%
      ggplot(aes(x=r2, fill= as.factor(pcr_plate))) +
      geom_histogram(binwidth=0.01, color='black', alpha=0.5) +
      labs(x= 'Sample R2', y= 'Frequency', fill= 'Plate', 
           title= 'Distribution of control barcode trend R2 values') + 
      theme_bw()
    
    pdf(file=paste(out, "control_barcode_trend_histogram.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
    print(trend_hist)
    dev.off()
  }
  
  ## Sample correlation -----
  print("generating sample_cor image")
  correlation_matrix= filtered_counts %>% dplyr::ungroup() %>% 
    dplyr::filter(!is.na(CCLE_name),
           (!trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type)) %>% 
    dplyr::mutate(log2_n = log2(n +1)) %>% 
    reshape2::dcast(CCLE_name~profile_id, value.var="log2_n") %>% 
    column_to_rownames("CCLE_name") %>% 
    cor(use="pairwise.complete.obs") 
  
  cp= correlation_matrix %>% reshape2::melt() %>% 
    ggplot() + geom_tile(aes(x=Var1, y=Var2, fill=value)) +
    labs(x="", y="", fill="correlation") +
    scale_fill_gradient(low="yellow", high="red") +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=5),
          axis.text.y = element_text(size=5))
  pdf(file=paste(out, "sample_cor.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
  print(cp)
  dev.off()
  
  ## Tech rep correlations ----
  # assumes that tech reps are the last component of profile_id
  if('tech_rep' %in% colnames(normalized_counts)) {
    tech_rep_cors= normalized_counts %>% reshape2::acast(DepMap_ID+Name~profile_id, value.var= 'log2_n') %>%
      cor(method='pearson', use='pair') %>% melt() %>%
      dplyr::mutate(bio_rep_id1= str_replace(Var1, ':\\d+$', ''),
                    tech_rep_val1= str_extract(Var1, '\\d+$') %>% as.numeric(),
                    bio_rep_id2= str_replace(Var2, ':\\d+$', ''),
                    tech_rep_val2= str_extract(Var2, '\\d+$') %>% as.numeric()) %>%
      dplyr::filter(bio_rep_id1 == bio_rep_id2, tech_rep_val2 > tech_rep_val1)
    
    tech_rep_hist= tech_rep_cors %>% ggplot(aes(x= value)) + geom_histogram() + 
      labs(x= 'Tech rep correlations', y= 'Frequency',
           title= 'Histogram of Tech rep corrleations') + 
      theme_bw()
    
    pdf(file=paste(out, "tech_rep_hist.pdf", sep="/"),
        width=sqrt(num_profiles), height=sqrt(num_profiles))
    print(tech_rep_hist)
    dev.off()
    
    if(max(unique(normalized_counts$tech_rep)) == 2) {
      print("generating tech rep correlations image")
   
      static_cols= c('project_code', 'CCLE_name', 'DepMap_ID', 'Name', 'cell_set')
      tech_reps_piv= normalized_counts %>% dplyr::mutate(bio_rep_id= str_replace(profile_id, ':\\d+$', '')) %>%
        dplyr::group_by_at(c(static_cols, 'bio_rep_id')) %>% dplyr::filter(n!=0, n()==2) %>% dplyr::ungroup() %>%
        pivot_wider(id_cols= all_of(c(static_cols, 'bio_rep_id')),
                    names_from= tech_rep, names_prefix= 'tech_rep', values_from= log2_n) %>%
        dplyr::group_by(bio_rep_id) %>%
        dplyr::mutate(r2= cor(tech_rep1, tech_rep2, use='p')^2,
                      type= ifelse(!is.na(CCLE_name), "cell line", "control barcode")) %>% dplyr::ungroup()
      
      tech_reps_plt= tech_reps_piv %>% dplyr::mutate(bio_rep_id= reorder(bio_rep_id, r2)) %>%
        ggplot(aes(x= tech_rep1, y= tech_rep2)) +
        geom_point(aes(color= type), alpha=0.75) +
        geom_smooth(method='lm', se=F, color='black', linewidth=0.5, linetype=2) +
        stat_correlation(mapping = use_label(c("R2", "n")))+ 
        facet_wrap(~bio_rep_id, scales= 'free') +
        labs(x="tech rep 1 log2(n)", y="tech rep 2 log2(n)") + theme_bw()
      
      pdf(file=paste(out, "tech_reps_plt.pdf", sep="/"),
          width=sqrt(num_profiles), height=sqrt(num_profiles))
      print(tech_reps_plt)
      dev.off()
    }
  }
  
  ## Bio rep correlations ----
  if('bio_rep' %in% colnames(normalized_counts)) {
    num_bio_reps= normalized_counts %>% 
      dplyr::filter((!trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type)) %>% 
      pull(bio_rep) %>% unique() %>% length()
    
    if(num_bio_reps > 1) {
      print("generating bio rep correlations image")
      # collapse tech reps taken from 'compute_l2fc'
      collapsed_tech_rep= normalized_counts %>%
        dplyr::filter(!(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type), !is.na(CCLE_name)) %>%
        dplyr::group_by_at(setdiff(names(.), c('pcr_plate','pcr_well', 'Name', 'log2_dose', 'cb_intercept',
                                               'profile_id', 'tech_rep', 'n', 'log2_n', 'normalized_n', 
                                               'log2_normalized_n', 'flag', count_col_name))) %>% 
        dplyr::summarise(mean_normalized_n = mean(!! rlang::sym(count_col_name)), 
                         num_tech_reps= n()) %>% dplyr::ungroup()
      collapsed_tech_rep$sig_id= do.call(paste,c(collapsed_tech_rep[sig_cols], sep=':'))
      
      bio_corr= collapsed_tech_rep %>% ungroup() %>% 
        filter(!is.na(CCLE_name),
               (!trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type)) %>% 
        mutate(plt_id= paste(sig_id, bio_rep, sep=':')) %>% 
        dcast(CCLE_name~plt_id, value.var="mean_normalized_n") %>% 
        column_to_rownames("CCLE_name") %>% 
        cor(use="pairwise.complete.obs")
      
      bio_corr_hm= bio_corr %>% melt() %>% ggplot() +
        geom_tile(aes(x=Var1, y=Var2, fill=value)) +
        labs(x="", y="", fill="correlation\nnorm_n") +
        scale_fill_gradientn(breaks= c(0, 0.5, 1), 
                             colours= c('white', 'yellow', 'red'),
                             limits=c(0,1), oob=squish) +
        theme(axis.text.x = element_text(angle=70, hjust=1, size=5),
              axis.text.y = element_text(size=5))
      
      pdf(file=paste(out, "bio_corr_hm.pdf", sep="/"),
          width=sqrt(num_profiles), height=sqrt(num_profiles))
      print(bio_corr_hm)
      dev.off()
      
      # bio cor histogram
      bio_rep_cor_long= bio_corr %>% melt() %>%
        dplyr::mutate(id1= str_replace(Var1, ':\\d+$', ''),
                      br1= str_extract(Var1, '\\d+$') %>% as.numeric(),
                      id2= str_replace(Var2, ':\\d+$', ''),
                      br2= str_extract(Var2, '\\d+$') %>% as.numeric()) %>%
        dplyr::filter(id1==id2, br2>br1)
      
      bio_rep_hist= bio_rep_cor_long %>% ggplot(aes(x= value)) + geom_histogram() + 
        labs(x= 'Bio rep correlations', y= 'Frequency',
             title= 'Histogram of bio rep corrleations') + 
        theme_bw()
      
      pdf(file=paste(out, "bio_rep_hist.pdf", sep="/"),
          width=sqrt(num_profiles), height=sqrt(num_profiles))
      print(bio_rep_hist)
      dev.off()
    }
  }
  
  # End _________________________ ----
  print('QC finishing')
}
