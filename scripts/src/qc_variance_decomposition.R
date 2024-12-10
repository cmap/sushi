require(data.table)
require(tidyverse)
require(magrittr) 
# The function currently outputs one dataframe that is in principle got values that are
# aggregated at different levels:
# 1. cell line/ CBC level withing a pool: var_log_fline.SEQ
# 2. cell line/CBC within a cell set: var_log_cl_in_pool.SEQ
# 3. cell pool/CBC pool within a cell set: var_log_fpool.SEQ
# these are all aggregated into one df that is returned by the function, but you can split them out if required.

### specific to a pool level or a cell line for each pcr plate

compute_variance_decomposition <- function(df, metric = 'log2_normalized_n',
                               negcon = "ctl_vehicle")
{
    
    ## the columns here don't always match ID_cols, SIG_cols etc so they are currently hardcoded.
    ## In principle each one could be passed as a different parameter?
    df = df %>% dplyr::mutate(pool_id=ifelse(!is.na(cb_name), "CTLBC", pool_id)) 
    ## make control barcodes have a pool id of CTLBC so that we can calculate properties for control barcodes as well
    
    var_log_fline.SEQ =df %>%
        filter(pert_type==negcon) %>%
        dplyr::group_by(cell_set, pcr_plate, pcr_well) %>%
        dplyr::mutate(tot_counts=sum(n+1)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(cell_set, pcr_plate, pcr_well,
                        depmap_id,cb_name) %>%
        dplyr::mutate(fcell_line=sum(n+1)/tot_counts) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(cell_set, pcr_plate,  
                        depmap_id, cb_name) %>%
        dplyr::summarise(var_log2_fline=var(log2(fcell_line)),
                         median_log2_fline=median(log2(fcell_line)),
                         mad_log2_fline=mad(log2(fcell_line)),
                         mean_log2_fline=mean(log2(fcell_line))) %>% 
        dplyr::ungroup()
    
    
    var_log_cl_in_pool.SEQ =  df %>%
        filter(pert_type==negcon) %>%
        dplyr::group_by( cell_set, pcr_plate, pcr_well,  
                        pool_id) %>%
        dplyr::mutate(tot_counts=sum(n+1)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by( cell_set, pcr_plate, pcr_well, 
                        depmap_id,cb_name, pool_id) %>%
        dplyr::mutate(fcl_in_pool=sum(n+1)/tot_counts) %>%
        dplyr::ungroup() %>%
        dplyr::group_by( cell_set, pcr_plate,  ,
                        depmap_id,cb_name, pool_id) %>%
        dplyr::summarise(var_log2_fcl_in_pool=var(log2(fcl_in_pool)),
                         mad_log2_fcl_in_pool=mad(log2(fcl_in_pool)),
                         median_log2_fcl_in_pool=median(log2(fcl_in_pool)),
                         mean_log2_fcl_in_pool=mean(log2(fcl_in_pool))) %>% 
        dplyr::ungroup()
    
    
    pwise_negcon_stats.SEQ = df %>%
        filter(pert_type==negcon) %>%
        dplyr::group_by( cell_set, pcr_plate, pcr_well,  ) %>%
        dplyr::mutate(tot_counts=sum(n+1)) %>% 
        dplyr::ungroup() %>%
        dplyr::group_by( cell_set, pcr_plate, pcr_well,  
                        pool_id) %>%
        dplyr::summarise(
            frac_reads=sum(n+1)/tot_counts) %>% 
        dplyr::ungroup() %>%
        dplyr::distinct()
    
    
    var_log_fpool.SEQ = pwise_negcon_stats.SEQ %>%
        dplyr::group_by( cell_set, pcr_plate,  
                        pool_id) %>% 
        dplyr::summarise(var_log2_frac_pool_reads=var(log2(frac_reads)),
                         mad_log2_frac_pool_reads=mad(log2(frac_reads)),
                         median_log2_frac_pool_reads=median(log2(frac_reads)),
                         mean_log2_frac_pool_reads=mean(log2(frac_reads))) %>% 
        dplyr::ungroup()
    
    
    var_decomp.SEQ <- dplyr::left_join(var_log_cl_in_pool.SEQ, var_log_fpool.SEQ,
                                       by=c("pool_id","cell_set", "pcr_plate")) %>%
        dplyr::left_join(var_log_fline.SEQ)

return(var_decomp.SEQ)
}


var_decomp.SEQ <- compute_variance_decomposition(normcts,
                                                 metric = 'log2_normalized_n', negcon = "ctl_vehicle")



