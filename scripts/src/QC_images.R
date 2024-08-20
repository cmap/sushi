#' validate_columns_exist
#' 
#' This function checks that a list of columns are present in a dataframe.
#' Columns that were not found in the dataframe are printed out.
#' 
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @returns Boolean
validate_columns_exist= function(selected_cols, df) {
  # Check that all of selected_columns are in df - base::setdiff(A, B) = A[!A %in% B].
  unmatched_cols= base::setdiff(selected_cols, colnames(df))
  
  if(length(unmatched_cols) > 0) {
    print('The following columns are missing: ')
    print(unmatched_cols)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Calculate index summaries
#' 
#' Generates some simple summaries for each unique index.
#' 
#' @param df A dataframe which must contain the column "n" which represents the count of a read.
#' @param index_col The name of the column contain the index barcodes as a string. This column must be present in "df".
#' @param valid_indices. A vector of all the valid indices for "index_col".
#' @returns A dataframe with the follow columns:
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

#' Calculate purity metrics
#' 
#' Create the qc table with index purity and cell line purity.
#' 
#' @param raw_counts_uncollapsed Dataframe output from nori.
#' @param raw_counts Raw counts dataframe outputed from collate_fastq_reads.
#' @param filtered_counts Filtered counts dataframe outputed from filter_raw_reads.
#' @param value_col String name of the counts column present all three dataframes.
#' @param file_path Location to write out the output.
#' @returns Writes out a QC_table to the file_path.
create_qc_table= function(raw_counts_uncollapsed, raw_counts, filtered_counts, value_col= 'n', file_path) {
  # Validation: Check that value_col is present in the three files.
  if(!validate_columns_exist(value_col, raw_counts_uncollapsed)) {
    stop(paste0('The column ', value_col, " was not detected in uncollapsed raw counts."))
  }
  if(!validate_columns_exist(value_col, raw_counts)) {
    stop(paste0('The column ', value_col, " was not detected in raw counts."))
  }
  if(!validate_columns_exist(value_col, filtered_counts)) {
    stop(paste0('The column ', value_col, " was not detected in filtered counts."))
  }
  
  # Calculate purities
  index_purity= sum(raw_counts[[value_col]]) / sum(raw_counts_uncollapsed[[value_col]])
  print(paste0('Index purity: ', round(index_purity, 4)))
  cell_line_purity= sum(filtered_counts[[value_col]]) / sum(raw_counts[[value_col]])
  print(paste0('Cell line purity: ', round(cell_line_purity, 4)))
  qc_table= data.frame(index_purity= index_purity, cell_line_purity= cell_line_purity)
  
  # Write out table
  print(paste0('Writing QC table out to ', file_path))
  qc_table %>% write.csv(file_path, row.names= FALSE, quote= FALSE)
}

#' Total counts barplot
#' 
#' Creates the total counts barplot with bars colored by the barcode type,
#' either a cell line barcode or control barcode.
#'
#' @param filtered_counts Filtered counts dataframe.
#' @param id_cols Vector of columns names that identify each sample.
#' @param facet_col String name of the column in filtered_counts to facet the plot.
#'                  This can be left as NA if there isn't a column to facet on.
#' @returns Returns a ggplot object.
create_total_counts_barplot= function(filtered_counts, id_cols, facet_col= NA) {
  # Validation: Check that id_cols and facet_col exist in filtered counts.
  if(!validate_columns_exist(na.omit(c(id_cols, facet_col)), filtered_counts)) {
    stop('Some input columns were not detected in filtered counts.')
  }
  
  # Sum up reads 
  total_counts= filtered_counts %>%
    dplyr::mutate(barcode_type= case_when(!is.na(CCLE_name) ~ 'cell line',
                                          !is.na(Name) ~ 'ctrl barcode')) %>%
    tidyr::unite(all_of(id_cols), col= 'sample_id', sep= ':', remove= FALSE, na.rm= FALSE) %>%
    dplyr::group_by(pick(all_of(na.omit(c('sample_id', facet_col, 'barcode_type'))))) %>%
    dplyr::summarise(total_counts= sum(n)) %>% dplyr::ungroup()
  
  total_counts_plot= total_counts %>% 
    ggplot(aes(x=sample_id, y=total_counts, fill=barcode_type)) +
    geom_col(alpha=0.75, position='identity') +
    geom_hline(yintercept= 10^4, linetype=2) + 
    {if(!is.na(facet_col)) facet_wrap(~.data[[facet_col]], scale= 'free_x')} +
    labs(x= "Sample constructed using id_cols", y="Total counts", fill= 'Barcode\ntype', 
         title= 'Filtered counts - unstacked') + 
    theme_bw() + theme(axis.text.x = element_text(angle=70, hjust=1))
  
  return(total_counts_plot)
}


#' Cell line recover barplot
#' 
#' Creates barplots of the cell lines recovered. The parameter "plot_type" can be used to plot the percentage or
#' the total cell line counts on teh y axis. The parameter "include_ctrl_bcs" can be used to include the control
#' barcodes in the cell line count. 
#' 
#' @param filtered_counts Filtered counts dataframe.
#' @param id_cols Vector of column names that identify each sample.
#' @param facet_col String name of the column in filtered_counts to facet the plot.
#' @param value_col String name of the column in filtered_counts that contains the counts.
#' @param counts_threshold Threshold used to determine low counts.
#' @param plot_type String of either "percent" or "count" to adjust the y axis to be either the percentage or the 
#'                  total number of cell lines.
#' @param include_ctrl_bcs Boolean. Set to TRUE if control barcodes are to be counted. 
#' @returns Returns a ggplot plot.
create_recovery_barplot= function(filtered_counts, id_cols, facet_col= NA, value_col= 'n', count_threshold, 
                                  plot_type= 'percent', include_ctrl_bcs= FALSE) {
  # Validation: Check that id_cols, facet_col, or value_col exist in filtered counts.
  if(!validate_columns_exist(na.omit(c(id_cols, facet_col, value_col)), filtered_counts)) {
    stop('Some input columns were not detected in filtered counts.')
  }
  
  # Filter out control barcodes if it is specified.
  if(include_ctrl_bcs == FALSE) {
    filtered_counts= filtered_counts %>% dplyr::filter(is.na(Name))
  }
  
  # Count number of cell lines/ barcodes for a detection group.
  recovery= filtered_counts %>%
    dplyr::add_count(pick(all_of(id_cols)), name= 'total_num_cls') %>%
    dplyr::mutate(detect_type= case_when(.data[[value_col]] == 0 ~ 'Not detected',
                                         .data[[value_col]] <= count_threshold ~ 'Low counts',
                                         .data[[value_col]] > count_threshold ~ 'Detected')) %>%
    dplyr::count(pick(all_of(c(id_cols, facet_col, 'detect_type', 'total_num_cls'))), name= 'num_cls_by_type') %>%
    tidyr::unite(all_of(id_cols), col= 'sample_id', sep= ':', remove= FALSE, na.rm= FALSE) %>%
    dplyr::mutate(percent= (num_cls_by_type / total_num_cls) * 100)
  
  # Set the y axis depending on the plot type.
  if(plot_type == 'count') {
    y_col= 'num_cls_by_type'
    y_text= 'Number of cell lines'
  } else {
    if(plot_type != 'percent') {
      print(paste0('Warning: ', plot_type, ' is not a valid plot type. Please use either count or percent.'))
      print('Defaulting to percent plot.')
    }
    y_col= 'percent'
    y_text= 'Percentage of cell lines recovered (%)'
  }
  
  # Create recovery plot.
  recov_plot= recovery %>%
    ggplot(aes(x= sample_id, y= .data[[y_col]], fill= reorder(detect_type, dplyr::desc(detect_type)))) +
    geom_col(alpha=0.75, position='stack') +
    {if(!is.na(facet_col)) facet_wrap(~.data[[facet_col]], scale= 'free_x')} +
    labs(x= "Sample constructed using id_cols", y= y_text, fill= '', title= 'Cell line recovery') + 
    theme_bw() + theme(axis.text.x = element_text(angle=70, hjust=1))
  
  return(recov_plot)
}

#' Control barcode scatter plot
#' 
#' Creates a scatter plot of the control barcodes.
#' 
#' @param normalized_counts Dataframe output from the normalize module.
#' @param id_cols Vector of column names that identify every PCR well.
#' @param value_col Name of the column that contains the values.
#' @returns Returns a ggplot object.
create_ctrlBC_scatterplots= function(normalized_counts, id_cols, value_col= 'log2_n') {
  # Validation: Check that id_cols and value_col exist in filtered counts.
  if(!validate_columns_exist(c(id_cols, value_col), normalized_counts)) {
    stop('Some input columns were not detected in normalized counts.')
  }
  
  # Detect norm_r2 and norm_mae. If columns do not exist, then roughly calculate those columns.
  if(any(!c('norm_r2', 'norm_mae') %in% colnames(normalized_counts))) {
    print('WARNING: Columns "norm_r2" and/or "norm_mae" were not detected in normalized_counts.', quote= FALSE)
    print('Calculating both columns - this method may not be as robust as the normalize module.')
    
    normalized_counts= normalized_counts %>% 
      dplyr::filter(!is.na(Name), control_barcodes %in% c("Y", "T", T), n != 0) %>%
      dplyr::group_by(pick(all_of(id_cols))) %>%
      dplyr::mutate(mean_y= mean(log2_dose),
                    residual2= (log2_dose - log2_normalized_n)^2,
                    squares2= (log2_dose - mean_y)^2,
                    norm_r2= 1 - sum(residual2) / sum(squares2),
                    norm_mae= median(abs(log2_dose- log2_normalized_n))) %>% ungroup()
  } 
  
  # Filter for just the control barcodes, create a profile_id for faceting, 
  # and determine the x and y positions for the r2 + mae label.
  cb_trend= normalized_counts %>% dplyr::filter(!is.na(Name), control_barcodes %in% c("Y", "T", T)) %>%
    tidyr::unite(all_of(id_cols), col= 'profile_id', sep= ':', remove= TRUE) %>%
    dplyr::group_by(profile_id) %>% dplyr::mutate(label_x_pos= min(.data[[value_col]]),
                                                  label_y_pos= max(log2_dose)) %>% dplyr::ungroup()
  
  # Create control barcode trend plot
  trend_scatter_plot= cb_trend %>% ggplot(aes(x= .data[[value_col]], y= log2_dose)) + 
    geom_point() +
    geom_abline(aes(slope=1, intercept= cb_intercept) , color='blue', alpha= 0.5) +
    geom_text(aes(x= label_x_pos, y= label_y_pos,
                  label= paste0('r2= ', round(norm_r2, 4), '\nmae= ', round(norm_mae, 4))), 
              hjust='inward', vjust='inward', alpha= 0.5) +
    facet_wrap(~profile_id, scales='free_x') +
    labs(title= 'Linear fit of control barcodes across all samples') + theme_bw()
  
  return(trend_scatter_plot)
}

#' Heatmap of correlations
#' 
#' Creates a correlation heatmap. A matrix of values is created from the input_df. The row_id_cols
#' are used identify each row and the col_id_cols are used to identify each column. The value_col is 
#' used to fill the matrix. Correlations are then computed.
#' 
#' @import tidyverse
#' @import WGCNA
#' @import reshape2
#' @param input_df Dataframe.
#' @param row_id_cols Vector of column names from input_df that identifies the cell lines. For example,
#'                    this can be "DepMap_ID", "CCLE_name" if only cell lines exist. It can also be 
#'                    "DepMap_ID", "CCLE_name", "Name" if control barcodes are also present.
#' @param col_id_cols Vector of column names from input_df that identifies the PCR wells or conditions.
#'                    For example, this can be "pcr_plate", "pcr_well" or a list of conditions like those in sig_cols.
#' @param value_col String name of the column in input_df to be used as the values.
#' @param cor_method WGCNA correlation method. This defaults to "pearson".
#' @returns Returns a ggplot object.
create_cor_heatmap= function(input_df, row_id_cols, col_id_cols, value_col,
                             cor_method= 'pearson') {
  
  # Validate that specified columns are in the dataframe.
  if(!validate_columns_exist(c(row_id_cols, col_id_cols, value_col), input_df)) {
    stop('Not all columns were detected in the input dataframe.')
  }
  
  # Create row and column names for pivoting to a matrix
  correlation_mx= input_df %>%
    tidyr::unite(all_of(row_id_cols), col= 'row_id', sep= ':', remove= TRUE) %>%
    tidyr::unite(all_of(col_id_cols), col= 'col_id', sep= ':', remove= TRUE)
  
  # Check that the row and column ids specify one value.
  validate_ids= correlation_mx %>% dplyr::group_by(row_id, col_id) %>%
    dplyr::filter(dplyr::n() > 1) %>% dplyr::ungroup()
  if(nrow(validate_ids) != 0) {
    print('The provide columns specify more than one value.')
    print(head(validate_ids))
    stop('Multiple values detected for a unique combination of "row_id_cols" and "col_id_cols".')
  }
  
  # Pivot and calculate correlations
  correlation_mx= correlation_mx %>% reshape2::acast(row_id~col_id, value.var= value_col) %>%
    WGCNA::cor(use= 'pairwise.complete.obs', method= cor_method)
  
  # Create heatmap 
  cor_heatmap= correlation_mx %>% reshape2::melt() %>% 
    ggplot(aes(x= Var1, y= Var2, fill= value)) + 
    geom_tile() +
    labs(x= '', y= '', fill= '', title= paste0('Correlations using ', value_col)) +
    scale_fill_gradientn(breaks= c(0, 0.5, 1), 
                         colours= c('blue', 'white','red'),
                         limits=c(0, 1), oob= scales::squish) +
    theme(axis.text.x = element_text(angle=70, hjust=1))
  
  return(cor_heatmap)
}

#' Scatter plots of two replicates
#' 
#' From a long table, creates scatter plots to two replicates.
#' 
#' @param input_df Dataframe.
#' @param cell_line_cols List of column names used to identify each cell line or control barcode.
#' @param replicate_group_cols List of column names that describe a group of similar conditions.
#' @param replicate_col Name of the column that specifies the replicate. This column should not be 
#'                      in replicate_group_cols!
#' @param value_col Name of the column in input_df that contains the values.
#' @param x_axis_rep String of the replicate identifier that should be on the x axis of the plot.
#' @param y_axis_rep String of the replicate identifier that should be on the y axis of the plot.
#' @returns Returns a ggplot object or NULL if all entries are filtered out.
create_replicate_scatterplots= function(input_df, cell_line_cols, replicate_group_cols, replicate_col, value_col,
                                        x_axis_rep= '1', y_axis_rep= '2') {
  # Validation: Check that input columns are present in the dataframe.
  if(!validate_columns_exist(c(cell_line_cols, replicate_group_cols, replicate_col, value_col), input_df)) {
    stop('Some input columns were not detected in normalized counts.')
  }
  
  # Validation: Check that replicate_col is not in replicate_group_cols.
  if(replicate_col %in% replicate_group_cols) {
    stop(paste0(replicate_col, ' should not be included in replicate_group_cols!'))
  }
  
  reps_piv= input_df %>% 
    tidyr::unite(all_of(replicate_group_cols), col= 'replicate_group', sep= ':', remove= TRUE, na.rm= FALSE) %>%
    dplyr::group_by(pick(all_of(c(cell_line_cols, 'replicate_group')))) %>%
    dplyr::filter(!is.na(.data[[replicate_col]]), .data[[replicate_col]] != '', .data[[value_col]] != 0,
                  dplyr::n() >= 2) %>% 
    dplyr::ungroup() 
  
  # Retun a null object if no entries pass the filter.
  if(nrow(reps_piv) == 0) {return(NULL)}
  
  reps_piv= reps_piv %>%
    pivot_wider(id_cols= all_of(c(cell_line_cols, 'replicate_group')),
                names_from= replicate_col, names_prefix= replicate_col, values_from= value_col) %>%
    dplyr::mutate(type= ifelse(!is.na(CCLE_name), "cell line", "control barcode")) %>% dplyr::ungroup()
  
  # Create names of the columns to plot on xy axes
  x_col_name= paste0(replicate_col, x_axis_rep)
  y_col_name= paste0(replicate_col, y_axis_rep)
  
  reps_scatter= reps_piv %>% dplyr::filter(!is.na(.data[[x_col_name]]), !is.na(.data[[y_col_name]])) %>%
    ggplot(aes(x= .data[[x_col_name]], y= .data[[y_col_name]])) +
    geom_abline(color='black', linewidth=0.5, linetype=2) +
    geom_point(aes(color= type), alpha=0.75) +
    ggpmisc::stat_correlation(mapping = use_label(c("R2", "n")))+ 
    facet_wrap(~replicate_group, scales= 'free') +
    labs(x= paste0(replicate_col, '1 ', value_col), y= paste0(replicate_col, '2 ', value_col),
         title= paste0('Scatter plot of ', replicate_col, ' with ', value_col)) +
    theme_bw()
  
  return(reps_scatter)
}

#' QC_images
#'
#' Takes in various pipeline outputs and generates 11 QC files.
#'
#' @param raw_counts_uncollapsed Dataframe output from nori. This is used to generate purity metrics and
#'                               the index summaries.
#' @param raw_counts Raw counts dataframe from the collate_fastq_reads modules. This is used to generate puritu metrics.
#' @param annotated_counts Annotated counts dataframe from the filter_raw_reads module.
#' @param filtered_counts Filtered counts dataframe from the filter_raw_reads module.
#' @param normalized_counts Normalized counts dataframe from the normalize module. This is an optional parameter.
#' @param l2fc L2FC dataframe from the compute_l2fc module. This is used for the bio_reps plot. 
#' @param sample_meta Dataframe of the sample metadata for the sequencing run.
#' @param CB_meta Dataframe of the control barcode metadata. This is only used for the CDF plot.
#' @param cell_set_meta Dataframe of the cell set metadata. This is only used for the CDF plot.
#' @param cell_line_cols Vector of sample meta column names used to describe a cell line or barcode.
#' @param id_cols Vector of sample meta column names used to identify each PCR well. 
#'                This defaults to "pcr_plate", "pcr_well".
#' @param sig_cols Vector of sample meta column names used to identify a unique treatment condition.
#' @param control_type String of how the negative controls are designated in the trt_type column in the sample_meta.
#' @param count_threshold Threshold for low read counts.
#' @param reverse_index2 Boolean set to TRUE if the sequencing involved the reverse complement workflow.
#' @param out Path to the directory to save the QC images.
#' @returns NA. QC images are written out to the specified folder.
QC_images= function(raw_counts_uncollapsed, raw_counts, 
                    annotated_counts, filtered_counts, normalized_counts= NA, l2fc, 
                    sample_meta, CB_meta, cell_set_meta,
                    cell_line_cols, 
                    id_cols= c('pcr_plate', 'pcr_well'), sig_cols,
                    control_type= 'negcon', count_threshold= 40, 
                    reverse_index2= FALSE, out = NA) {
  require(tidyverse)
  require(magrittr)
  require(reshape2)
  require(scales)
  
  if(is.na(out)) {
    out = getwd()
  }
  
  # Some preprocessing ----
  skipped_qcs= c() # empty vector to collect potential errors
  num_profiles = annotated_counts %>% dplyr::distinct(pick(all_of(id_cols))) %>% nrow()
  
  # Detect control barcodes
  cb_check= sample_meta %>%
    dplyr::filter(control_barcodes %in% c("Y", "T", T),
                  !(trt_type %in% c("empty", "", "CB_only")) & !is.na(trt_type))
  contains_cbs= ifelse(nrow(cb_check)!= 0, T, F)
  
  # Count number of cell lines in each cell set
  num_cls_in_set= cell_set_meta %>% dplyr::filter(cell_set %in% unique(sample_meta$cell_set))
  # Add cell_cets that are 'missing' or are strings of LUAs
  if(nrow(num_cls_in_set) != length(unique(sample_meta$cell_set))) {
    sets_not_in_meta= sample_meta %>% dplyr::filter(!cell_set %in% cell_set_meta$cell_set) %>%
      dplyr::pull(cell_set) %>% unique() %>% sort()
    sets_to_add_df= data_frame(cell_set= sets_not_in_meta, members= sets_not_in_meta)
    num_cls_in_set= dplyr::bind_rows(num_cls_in_set, sets_to_add_df)
  }
  num_cls_in_set %<>% dplyr::mutate(expected_num_cl= str_split(members, ';')) %>%
    tidyr::unnest(cols= expected_num_cl) %>% dplyr::group_by(cell_set) %>% 
    dplyr::summarize(expected_num_cl= length(unique(expected_num_cl))) %>% dplyr::ungroup()
  #
  
  # Sequencing QCs ____________________ ----
  ## Purity metrics ----
  # call this function
  print('1. Generating QC table ...')
  create_qc_table(raw_counts_uncollapsed, raw_counts, filtered_counts,
                  value_col= 'n', file_path= paste0(out, '/QC_table.csv'))
  #
  
  ## Index count summaries ----
  print("2. Generating index counts tables ...")
  # Check that "IndexBarcode1" and "index_1" columns are present.
  # If so, calculate index summary and write out.
  if('index_1' %in% colnames(sample_meta) & 'index_1' %in% colnames(raw_counts_uncollapsed)) {
    expected_index1= unique(sample_meta$index_1)
    index1_counts= get_index_summary(raw_counts_uncollapsed, 'index_1', expected_index1)
    index1_counts %>% write.csv(file= paste(out, 'index1_counts.csv', sep='/'), row.names=F)
  } else {
    print('Column "index_1" not detected. Skipping index 1 summaries ...', quote= FALSE)
  }
  
  # Do the same for index 2.
  # Reverse index 2 barcodes if needed.
  if(reverse_index2) {
    print("Reverse-complementing index 2 barcode.")
    sample_meta$index_2= chartr("ATGC", "TACG", stringi::stri_reverse(sample_meta$index_2))
  }
  
  if('index_2' %in% colnames(sample_meta) & 'index_2' %in% colnames(raw_counts_uncollapsed)) {
    expected_index2= unique(sample_meta$index_2)
    index2_counts= get_index_summary(raw_counts_uncollapsed, 'index_2', expected_index2)
    index2_counts %>% write.csv(file= paste(out, 'index2_counts.csv', sep='/'), row.names=F)
  } else {
    print('Column "index_2" not detected. Skipping index 2 summaries ...',  quote= FALSE)
  }
  
  ## Total counts ----
  print("3. Generating total_counts image ...")
  potential_error= base::tryCatch({
    tc= create_total_counts_barplot(filtered_counts, id_cols, facet_col= 'pcr_plate')
    
    pdf(file=paste(out, "total_counts.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
    print(tc)
    dev.off()
    rm(tc)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the total counts barplot. Skipping this output ...') 
    return('QC table')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  # Assay QCs _________________________ ----
  ## Cell lines recovered ----
  print("4. Generating cell_lines_present image ...")
  potential_error= base::tryCatch({
    cl_rec= create_recovery_barplot(filtered_counts, id_cols= id_cols, facet_col= 'pcr_plate', 
                                    count_threshold= count_threshold, plot_type= 'percent')
    
    pdf(file=paste(out, "cell_lines_present.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles))
    print(cl_rec)
    dev.off()
    rm(cl_rec)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the recovery barplot. Skipping this output ...')
    return('CL recovery')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## Contaminants ----
  print('5. Generating contaminate cell lines ...')
  potential_error= base::tryCatch({
    contams= annotated_counts %>% dplyr::filter(expected_read==F) %>%
      dplyr::mutate(barcode_id= ifelse(is.na(CCLE_name), Name, CCLE_name)) %>%
      dplyr::group_by(forward_read_cl_barcode, barcode_id) %>% 
      dplyr::summarise(num_wells= n(), median_n=median(n), max_n= max(n)) %>% ungroup() %>%
      dplyr::arrange(desc(num_wells))
    
    contams %>% write.csv(file= paste(out, 'contam_cell_lines.csv', sep='/'), row.names=F)
    rm(contams)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the contaminants file. Skipping this output ...') 
    return('contam_cell_lines.csv')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## Contaminates for ursula ----
  print('6. Generating contaminate reads for Ursula ...')
  potential_error= base::tryCatch({
    pcr_locations= c('pcr_plate', 'pcr_well')
    
    # Validation: Check that the PCR columns are present in raw_counts.
    if(!validate_columns_exist(pcr_locations, raw_counts)) {
      stop('pcr_plate and pcr_well are required raw_counts.csv for this to work.')
    }
      
    # count number of wells a cell_set appears in.
    pcr_plate_map= sample_meta %>% dplyr::distinct(pick(any_of(c(pcr_locations, 'cell_set')))) %>%
      dplyr::group_by(pcr_plate) %>% dplyr::mutate(num_wells_in_plate= dplyr::n()) %>% dplyr::ungroup() %>%
      dplyr::group_by(cell_set) %>% dplyr::mutate(num_wells_in_set= dplyr::n()) %>% dplyr::ungroup()

    # index filter and identify reads as mapped or not
    sequencing_filter= raw_counts %>%
      dplyr::mutate(mapped= forward_read_cl_barcode %in% unique(annotated_counts$forward_read_cl_barcode))

    # total counts per well - used to calculate fractions
    counts_per_well= sequencing_filter %>% dplyr::group_by(pick(all_of(pcr_locations))) %>%
      dplyr::summarise(well_total_n= sum(n)) %>% dplyr::ungroup()

    # mapped contaminates to bind
    mapped_contams= annotated_counts %>% dplyr::filter(!expected_read) %>%
      dplyr::mutate(barcode_name= ifelse(is.na(CCLE_name), Name, CCLE_name)) %>%
      dplyr::select(all_of(c(pcr_locations, 'forward_read_cl_barcode', 'n', 'barcode_name')))
    
    contam_reads= sequencing_filter %>% dplyr::filter(mapped == FALSE) %>% dplyr::select(-mapped) %>%
      dplyr::bind_rows(mapped_contams) %>%
      dplyr::left_join(counts_per_well, by= pcr_locations) %>%
      dplyr::left_join(pcr_plate_map, by= pcr_locations) %>%
      # filter out barcodes that only appear in one well
      dplyr::group_by(forward_read_cl_barcode) %>% dplyr::filter(dplyr::n() >1) %>% dplyr::ungroup() %>%
      # number of wells in a pcr plate a barcode is detected in
      dplyr::group_by(forward_read_cl_barcode, pcr_plate) %>%
      dplyr::mutate(num_wells_detected_plate= dplyr::n()) %>% dplyr::ungroup() %>%
      # number of wells in a cell set a barcode is detected in
      dplyr::group_by(forward_read_cl_barcode, cell_set) %>%
      dplyr::mutate(num_wells_detected_set= dplyr::n()) %>% dplyr::ungroup() %>%
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
    
    # write out
    contam_reads %>% write.csv(paste0(out, 'contam_reads.csv'), row.names=F)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the contams for UW file. Skipping this output ...')
    return('contam for UW')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## Cumulative counts by lines in negcons ----
  print("7. Generating cummulative image ...")
  potential_error= base::tryCatch({
    cdf= filtered_counts %>% dplyr::filter(trt_type == control_type) %>%
      dplyr::left_join(num_cls_in_set, by= "cell_set") %>%
      dplyr::mutate(expected_num_cl= ifelse(control_barcodes, expected_num_cl + length(unique(CB_meta$Name)),
                                            expected_num_cl)) %>% # add CBs to expected_num_cl if there are CBs
      tidyr::unite(all_of(id_cols), col= 'profile_id', sep= ':', remove= FALSE) %>%
      dplyr::group_by(pcr_plate, pcr_well, profile_id, expected_num_cl) %>% 
      dplyr::mutate(total_counts= sum(n), pct_counts= n/total_counts,) %>% dplyr::arrange(-n) %>% 
      dplyr::mutate(cum_pct= cumsum(pct_counts), rank= row_number(),
                    rank_pct= rank/expected_num_cl) %>% dplyr::ungroup()
    
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
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the cdf plot. Skipping this output ...') 
    return('cdf plot')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## Control barcode trends ----
  if(contains_cbs & is.data.frame(normalized_counts)) {
    print("8. Generating control_barcode_trend image")
    potential_error= base::tryCatch({
      trend_sc= create_ctrlBC_scatterplots(normalized_counts, id_cols, value_col= 'log2_n')
      
      pdf(file=paste(out, "control_barcode_trend.pdf", sep="/"),
          width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
      print(trend_sc)
      dev.off()
      rm(cb_trend, trend_sc)
    }, error= function(e) {
      print(e)
      print('Encountered an error when creating the CB trends plot. Skipping this output ...') 
      return('cb trend')
    })
    
    # Collect returned string if an error occurred
    if(!is.null(potential_error)) {
      skipped_qcs = c(skipped_qcs, potential_error)
    }
  } else {
    print('8. No control barcodes detected. Skipping control_barcode_trend image.')
  }
  
  ## Sample correlation -----
  print("9. Generating sample_cor image ...")
  potential_error= base::tryCatch({
    cor_df= filtered_counts %>% 
      dplyr::filter(!is.na(DepMap_ID), !is.na(trt_type), !trt_type %in% c("empty", "", "CB_only")) %>%
      dplyr::mutate(log2_n= log2(n + 1))
    cp= create_cor_heatmap(input_df= cor_df,
                           row_id_cols= c('DepMap_ID'),
                           col_id_cols= c(sig_cols, id_cols),
                           value_col= 'log2_n')
    
    pdf(file=paste(out, "sample_cor.pdf", sep="/"),
        width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
    print(cp)
    dev.off()
    rm(correlation_matrix, cp)
  }, error= function(e) {
    print(e)
    print('Encountered an error when creating the sample_cor figure. Skipping this output ...')
    return('sample_cor')
  })
  
  # Collect returned string if an error occurred
  if(!is.null(potential_error)) {
    skipped_qcs = c(skipped_qcs, potential_error)
  }
  
  ## Tech rep correlations ----
  if(is.data.frame(normalized_counts) & 'tech_rep' %in% colnames(normalized_counts)) {
    # Check if there are more at least two tech reps
    unique_tech_reps= na.omit(unique(normalized_counts$tech_rep))
    
    if(length(unique_tech_reps) >= 2) {
      print("10. Generating tech rep correlations image ...")
      # Set up replicate groups depending "bio_rep" column
      if('bio_rep' %in% colnames(normalized_counts) & !'bio_rep' %in% sig_cols) {
        replicate_group_cols= c(sig_cols, 'bio_rep')
      } else {
        replicate_group_cols= sig_cols
      }
      
      # Handle cases if control barcodes are used.
      if('Name' %in% colnames(normalized_counts)) {
        unique_cell_line_cols= c(cell_line_cols, 'Name')
      } else {
        unique_cell_line_cols= cell_line_cols
      }
      
      # Create replicate scatter plot
      potential_error= base::tryCatch({
        tech_reps_plt= create_replicate_scatterplots(input_df= normalized_counts, 
                                                     cell_line_cols= unique_cell_line_cols, 
                                                     replicate_group_cols= replicate_group_cols, 
                                                     replicate_col= 'tech_rep', 
                                                     value_col= 'log2_n')
        if(!is.null(tech_reps_plt)) {
          pdf(file=paste(out, "tech_reps_plt.pdf", sep="/"),
              width=sqrt(num_profiles), height=sqrt(num_profiles))
          print(tech_reps_plt)
          dev.off()
        } else {
          print('No technical replicates detected - skipping plot.')
        }
      }, error= function(e) {
        print(e)
        print('Encountered an error when creating the tech_reps_plt figure. Skipping this output ...')
        return('tech_reps_plt')
      })
      
      # Collect returned string if an error occurred
      if(!is.null(potential_error)) {
        skipped_qcs = c(skipped_qcs, potential_error)
      }
      
    } else {
      print('10. No technical replicates detected. Skipping tech_reps scatter plot.')
    }
  } else {
    print('10. No technical replicates detected. Skipping tech_reps scatter plot.')
  }
  
  ## Bio rep correlations ----
  if('bio_rep' %in% colnames(l2fc)) {
    unique_bio_reps= na.omit(unique(l2fc$bio_rep))
    
    if(length(unique_bio_reps) >= 2) {
      l2fc_with_log2= l2fc %>% dplyr::mutate(log2_mean_normalized_n= log2(mean_normalized_n))
      
      # Bio replicate scatter plots
      # bio_reps_plt= create_replicate_scatterplots(input_df= l2fc_with_log2s, 
      #                                           cell_line_cols= cell_line_cols, 
      #                                           replicate_group_cols= sig_cols, 
      #                                           replicate_col= 'bio_rep', 
      #                                           value_col= 'log2_mean_normalized_n')
      # if(!is.null(bio_reps_plt)) {
      #   pdf(file=paste(out, "bio_reps_plt.pdf", sep="/"),
      #       width=sqrt(num_profiles), height=sqrt(num_profiles))
      #   print(bio_reps_plt)
      #   dev.off()
      # } else {
      #   print('No technical replicates detected - skipping plot.')
      # }
      
      # Bio replicate heatmap
      print("11. Generating bio rep correlations heatmap ...")
      potential_error= base::tryCatch({
        bio_corr_hm= create_cor_heatmap(input_df= l2fc_with_log2, 
                                        row_id_cols= cell_line_cols, 
                                        col_id_cols= c(sig_cols, 'bio_rep'), 
                                        value_col= 'l2fc',
                                        cor_method= 'pearson') 
        pdf(file=paste(out, "bio_corr_hm.pdf", sep="/"),
            width=sqrt(num_profiles), height=sqrt(num_profiles))
        #print(bio_corr_hm)
        #dev.off()
      }, error= function(e) {
        print(e)
        print('Encountered an error when creating the bio_corr_hm figure. Skipping this output ...') 
        return('bio_corr_hm')
      })
      
      # Collect returned string if an error occurred
      if(!is.null(potential_error)) {
        skipped_qcs = c(skipped_qcs, potential_error)
      }
      
    } else {
      print('11. No biological replicates detected. Skipping bio_rep heatmap.')
    }
  }
  
  # End _________________________ ----
  print('QC finishing')
  if(length(na.omit(skipped_qcs)) != 0) {
    print(paste0('WARNING: The following ', length(skipped_qcs), ' QCs encountered errors and were skipped - '))
    print(na.omit(skipped_qcs))
  } else {
    print('No errors encountered.')
  }
}
