#' Calculate index summaries
#' 
#' Generates some simple summaries for each unique index.
#' 
#' @param df A dataframe which must contain the column "n" which represents the count of a read.
#' @param index_col The name of the column contain the index barcodes as a string. This column must be present in "df".
#' @param valid_indices. A vector of all the valid indices for "index_col".
#' @return A dataframe with the follow columns:
#'         - index_col: String, The column containing the index barcodes.
#'         - idx_n: Numeric, Number of reads associated with a specific index barcode.
#'         - fraction: Numeric, "idx_n" divided by the total number of reads in the run.
#'         - expected: Boolean, True if the index barcode is in "valid_indices" otherwise False.
#'         - contains_n: Boolean, True if the index barcode contains "N" in its sequence, otherwise False.
#'         - lv_dist: Numeric, Edit distance from a valid index barcode.
#'         - ham_dist: Numeric, Hamming distance from a valid index barcode.
get_index_summary= function(df, index_col, valid_indices) {
  output_summary= df %>% dplyr::group_by(pick(all_of(index_col))) %>% 
    dplyr::summarise(idx_n= sum(n)) %>% dplyr::ungroup() %>%
    dplyr::mutate(fraction= round(idx_n/sum(idx_n), 5),
                  expected= ifelse(.[[index_col]] %in% valid_indices, T, F),
                  contains_n= ifelse(grepl('N', .[[index_col]]), T, F),
                  lv_dist= apply(stringdist::stringdistmatrix(.[[index_col]], valid_indices, method="lv"), 
                                 1, min),
                  ham_dist= apply(stringdist::stringdistmatrix(.[[index_col]], valid_indices, method="hamming"), 
                                  1, min)) %>%
    dplyr::arrange(desc(fraction))
  return(output_summary)
}

#'  QC_images
#'
#'  Takes in the metadata, raw counts, annotated counts, and normalized counts to generate some QC images.
#'  
#' @param sample_meta - sample metadata
#' @param annotated_counts - dataframe of annotated readcounts that must include the following columns:
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
#' @param count_col_names - which counts to plot
#' @param control_type - how the negative controls are designated in the trt_type column in the sample metadata
#' @param count_threshold - threshold for low counts
#' @param reverse_index2 reverse index 2 if newer sequencers are used.
#' @return - NA, QC images are written out to the specified folder
#' @export

QC_images = function(sample_meta, raw_counts, annotated_counts, normalized_counts= NA,
                     CB_meta, cell_set_meta, out = NA, sig_cols, count_col_name= 'normalized_n',
                     control_type, count_threshold= 40, reverse_index2= FALSE) {
  require(tidyverse)
  require(magrittr)
  require(reshape2)
  require(scales)
  
  if(is.na(out)) {
    out = getwd()
  }
  num_profiles = annotated_counts$profile_id %>% unique() %>% length()
  
  # Some preprocessing ----
  # Reverse index 2 barcodes
  if(reverse_index2) {
    print("Reverse-complementing index 2 barcode.")
    sample_meta$index_2= chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$index_2))
  }
  
  # Detect control barcodes
  cb_check= sample_meta %>%
    dplyr::filter(control_barcodes %in% c("Y", "T", T),
                  !(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type))
  contains_cbs= ifelse(nrow(cb_check)!= 0, T, F)
  #
  
  # Sequencing QCs ____________________ ----
  ## Index count summaries ----
  print("Generating index counts tables")
  # Check that "IndexBarcode1" and "index_1" columns are present.
  # If so, calculate index summary and write out.
  if('index_1' %in% colnames(sample_meta) & 'index_1' %in% colnames(raw_counts)) {
    expected_index1= unique(sample_meta$index_1)
    index1_counts= get_index_summary(raw_counts, 'index_1', expected_index1)
    index1_counts %>% write.csv(file= paste(out, 'index1_counts.csv', sep='/'), row.names=F)
  } else {
    print('Column "index_1" not detected. Skipping index 1 summaries ...')
  }
  
  # Do the same for index 2.
  if('index_2' %in% colnames(sample_meta) & 'index_2' %in% colnames(raw_counts)) {
    expected_index2= unique(sample_meta$index_2)
    index2_counts= get_index_summary(raw_counts, 'index_2', expected_index2)
    index2_counts %>% write.csv(file= paste(out, 'index2_counts.csv', sep='/'), row.names=F)
  } else {
    print('Column "index_2" not detected. Skipping index 2 summaries ...')
  }
  
  ## Total counts ----
  print("generating total_counts image")
  total_counts= annotated_counts %>% dplyr::filter(expected_read) %>%
    mutate(type = ifelse(!is.na(CCLE_name), "cell line", "control barcode"),
           sample_id= paste(pcr_plate, pcr_well, sep='_')) %>%
    group_by(pcr_plate, pcr_well, sample_id, profile_id, type) %>% 
    dplyr::summarise(total_counts = sum(n))
  
  tc= total_counts %>% ggplot() +
    geom_col(aes(x=sample_id, y=total_counts, fill=type), alpha=0.75, position='identity') +
    geom_hline(yintercept= 10^4, linetype=2) + 
    facet_wrap(~pcr_plate, scale= 'free_x') +
    labs(x="PCR location", y="total counts", fill="", title= 'Raw counts - unstacked') + theme_bw() +
    theme(axis.text.x = element_text(angle=70, hjust=1, size=5)) 
  
  pdf(file=paste(out, "total_counts.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
  print(tc)
  dev.off()
  rm(total_counts, tc)
  
  # Assay QCs _________________________ ----
  ## Cell lines recovered ----
  print("generating cell_lines_present image")
  recovery= annotated_counts %>% dplyr::filter(expected_read, !is.na(CCLE_name)) %>%
    dplyr::select(!any_of(c('members'))) %>%
    dplyr::left_join(cell_set_meta, by="cell_set", relationship= 'many-to-one') %>%
    dplyr::mutate(members= if_else(is.na(members), cell_set, members), # for custom cell sets
                  expected_num_cl= as.character(members) %>% purrr::map(strsplit, ";") %>% 
                    purrr::map(`[[`, 1) %>% purrr::map(length) %>% as.numeric(),
                  count_type= ifelse(n > count_threshold, 'Detected', 'Low'),
                  count_type= ifelse(n==0, 'Missing', count_type)) %>%
    dplyr::count(profile_id, pcr_plate, pcr_well, count_type, expected_num_cl, name= 'count') %>%
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
  rm(recovery, cl_rec)
  
  ## Contaminants ----
  contams= annotated_counts %>% dplyr::filter(expected_read==F) %>%
    dplyr::mutate(barcode_id= ifelse(is.na(CCLE_name), Name, CCLE_name)) %>%
    dplyr::group_by(forward_read_cl_barcode, barcode_id) %>% 
    dplyr::summarise(num_wells= n(), median_n=median(n), max_n= max(n)) %>% ungroup() %>%
    dplyr::arrange(desc(num_wells))
  
  contams %>% write.csv(file= paste(out, 'contams.csv', sep='/'), row.names=F)
  rm(contams)
  
  # Work in progress - Updating contaminant output with new filtered counts
  run_below= F
  if(run_below) {
    # Determine which seq cols are present.
    # sample_meta sequencing columns
    sm_seq_cols= c('flowcell_name', 'flowcell_lane', 'IndexBarcode1', 'IndexBarcode2')
    # raw_counts sequencing columns
    rc_seq_cols= c('flowcell_name', 'flowcell_lane', 'index_1', 'index_2')
    
    mapped_reads= annotated_counts %>% 
      dplyr::distinct(pick(any_of(c(sm_seq_cols, 'forward_read_cl_barcode', 'expected')))) %>%
      dplyr::mutate(mapped=T)
    
    meta_joining_vector= list()
    for(item in seq_cols) {
      if(item=='IndexBarcode1') {
        meta_joining_vector[['index_1']]= item
      } else if(item=='IndexBarcode2') {
        meta_joining_vector[['index_2']]= item
      } else {
        meta_joining_vector[[item]]= item
      }
    }
    meta_joining_vector= unlist(meta_joining_vector)
    
    # unique combos of seq_cols
    unique_seq_col_vals= sample_meta %>% dplyr::distinct(pick(any_of(sm_seq_cols)))
    
    # get unmapped reads
    unmapped_read= raw_counts %>% 
      dplyr::semi_join(unique_seq_col_vals, by= meta_joining_vector) %>%
      dplyr::mutate(mapped= ifelse(forward_read_cl_barcode %in% c(cell_line_meta$Sequence, CB_meta$Sequence), T, F))
    
    dplyr::filter(index_1 %in% sample_meta$IndexBarcode1, index_2 %in% sample_meta$IndexBarcode2) %>%
      dplyr::left_join(mapped_reads, by= c('index_1', 'index_2', 'forward_read_cl_barcode')) %>%
      tidyr::replace_na(list(mapped= F)) %>% dplyr::filter(mapped== F) %>% dplyr::select(-mapped)
    
    # map of index barcodes to pcr_plate location
    pcr_plate_map= sample_meta %>%
      dplyr::distinct(pick(any_of(c('IndexBarcode1', 'IndexBarcode2', 'pcr_plate', 'pcr_well', 'cell_set')))) %>%
      dplyr::group_by(pcr_plate) %>% dplyr::mutate(num_wells_in_plate= dplyr::n()) %>% dplyr::ungroup() %>%
      dplyr::group_by(cell_set) %>% dplyr::mutate(num_wells_in_set= dplyr::n()) %>% dplyr::ungroup()
    
    # total counts per well - used to calculate fractions
    counts_per_well= raw_counts %>%
      dplyr::filter(index_1 %in% sample_meta$IndexBarcode1, index_2 %in% sample_meta$IndexBarcode2) %>%
      dplyr::group_by(index_1, index_2) %>% dplyr::summarise(well_total_n= sum(n)) %>% dplyr::ungroup()
    
    # mapped contaminates to bind
    mapped_contams= annotated_counts %>% dplyr::filter(!expected_read) %>%
      dplyr::mutate(barcode_name= ifelse(is.na(CCLE_name), Name, CCLE_name)) %>%
      dplyr::select(index_1, index_2, forward_read_cl_barcode, barcode_name, n)
    
    # put everything together
    contam_reads= unmapped_read %>% dplyr::bind_rows(mapped_contams) %>%
      dplyr::left_join(counts_per_well, by= c('index_1', 'index_2')) %>%
      dplyr::left_join(pcr_plate_map, by= join_by('index_1'=='IndexBarcode1', 'index_2'=='IndexBarcode2')) %>%
      # filter out barcodes that only appear in one well
      dplyr::group_by(forward_read_cl_barcode) %>% dplyr::filter(dplyr::n() >1) %>% dplyr::ungroup() %>%
      # number of wells in a pcr plate a barcode is detected in
      dplyr::group_by(forward_read_cl_barcode, pcr_plate) %>%
      dplyr::mutate(num_wells_detected_plate= n()) %>% dplyr::ungroup() %>%
      # number of wells in a cell set a barcode is detected in
      dplyr::group_by(forward_read_cl_barcode, cell_set) %>%
      dplyr::mutate(num_wells_detected_set= n()) %>% dplyr::ungroup() %>%
      # determine if contamination is project, plate, or set
      dplyr::group_by(forward_read_cl_barcode) %>%
      dplyr::mutate(num_wells_detected= dplyr::n(),
                    project_code= unique(sample_meta$project_code),
                    fraction= n/well_total_n,
                    type1= ifelse(sum(num_wells_detected== nrow(pcr_plate_map))>1, 'project_contam', NA),
                    type2= ifelse(sum(num_wells_detected== num_wells_detected_plate & 
                                        num_wells_detected_plate == num_wells_in_plate)>1, 'plate_contam', NA),
                    type3= ifelse(sum(num_wells_detected == num_wells_detected_set &
                                        num_wells_detected_set== num_wells_in_set)>1, 'set_contam', NA)) %>%
      dplyr::ungroup() %>%
      tidyr::unite(scope, all_of(c('type1', 'type2', 'type3')), sep=',', remove = T, na.rm = T) %>%
      dplyr::group_by(project_code, forward_read_cl_barcode, barcode_name, scope, num_wells_detected) %>%
      dplyr::summarise(min_n= min(n), med_n= median(n), max_n= max(n),
                       min_fraction= min(fraction), med_fraction= median(fraction), max_fraction=max(fraction)) %>%
      dplyr::arrange(desc(max_fraction))
    # GCCTATTAATCATTACTACTAATC
    contam_reads %>% write.csv(paste0(folder_path, 'contam_reads.csv'), row.names=F)
  }
  
  ## Cumulative counts by lines in negcons ----
  print("generating cummulative image")
  cdf= annotated_counts %>% dplyr::filter(expected_read) %>% 
    dplyr::filter(trt_type== control_type) %>% # filter for only negative controls
    dplyr::select(!any_of(c('members'))) %>%
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
  rm(cdf, mark50, mark95, cdf_plot)
  
  ## Control barcode trends ----
  if(contains_cbs & is.data.frame(normalized_counts)) {
    print("generating control_barcode_trend image")
    
    # calculate r2 and mae if columns do not exist
    if (!'r2' %in% colnames(normalized_counts) | !'r2' %in% colnames(normalized_counts)) {
      cb_trend= normalized_counts %>%
        dplyr::filter(control_barcodes %in% c("Y", "T", T),
                      !(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type),
                      !is.na(Name)) %>%
        dplyr::group_by(profile_id) %>%
        dplyr::mutate(mean_y= mean(log2_dose),
                      residual2= (log2_dose - log2_normalized_n)^2,
                      squares2= (log2_dose - mean_y)^2,
                      r2= 1 - sum(residual2)/sum(squares2),
                      mae= median(abs(log2_dose- log2_normalized_n))) %>% ungroup()
    } else {
      cb_trend= normalized_counts %>% dplyr::filter(!is.na(Name))
    }
    
    trend_sc= cb_trend %>% dplyr::mutate(profile_id= reorder(profile_id, desc(mae))) %>%
      ggplot(aes(x=log2_n, y=log2_dose)) + geom_point() +
      geom_abline(aes(slope=1, intercept= cb_intercept) , color='blue') +
      geom_text(aes(x= min(log2_n), y= desc(sort(unique(log2_dose)))[1], 
                    label= paste('r2=', round(r2, 4), '\nmae=', round(mae, 4), sep='')), 
                hjust='inward', vjust='inward') +
      facet_wrap(~profile_id, scales= 'free_x') +
      labs(x= 'log2(n)', y= 'log2(dose)') + theme_bw()
    
    pdf(file=paste(out, "control_barcode_trend.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
    print(trend_sc)
    dev.off()
    rm(cb_trend, trend_sc)
  }
  
  ## Sample correlation -----
  print("generating sample_cor image")
  correlation_matrix= annotated_counts %>% dplyr::filter(expected_read) %>% 
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
  rm(correlation_matrix, cp)
  
  ## Tech rep correlations ----
  # assumes that tech reps are the last component of profile_id
  if('tech_rep' %in% colnames(normalized_counts)) {
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
        reshape2::acast(CCLE_name~plt_id, value.var="mean_normalized_n") %>% 
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
    }
  }
  
  # End _________________________ ----
  print('QC finishing')
}
