options(cli.unicode = FALSE)
# Generates the QC table for EPS.
# It reads in normalized counts.csv  and uses  count_threshold and days_to_drop parameters to:
# flag cell lines  that do not meet the count_threshold in negcon wells on the day that isn't dropped (only a single day can be presented on the portal)

##Assay pool information will need to be added when merging has been implemented in pipeline
## so now QC table has cell set information rather than assay pool information
# import necessary libraries and functions
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(data.table))

#---- Read arguments ----
# initialize parser
parser <- ArgumentParser()
# specify our desired options
parser$add_argument("-b", "--build_dir", default="", help="Input Directory with data from whole screen")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")
parser$add_argument("-n", "--name", default="", help = "Build name. Default is none")
parser$add_argument("--control_type", default = "ctl_vehicle", help = "how negative control wells are distinguished in the trt_type column")
parser$add_argument("--count_threshold", default= 40, help = "Low counts threshold")
parser$add_argument("-d","--days",default="", help='Day timepoints to drop from output data separated by commas.')

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

base_dir <- args$build_dir
out_dir <- args$out
build_name <- args$name
control_type <- args$control_type
days_to_drop.arg <- args$days
count_threshold <- as.numeric(args$count_threshold) #default is 40 
days_to_drop <- as.numeric(unlist(strsplit(days_to_drop.arg, ","))) ## convert from string to numeric vector

if (!dir.exists(out_dir)) {dir.create(out_dir, recursive = T)} ## create output dir if it doesn't exist

## read data
norm_count_files <- list.files(base_dir, pattern = "normalized_counts.csv", full.names = T)
if (length(norm_count_files) == 1) {
    norm_counts <- data.table::fread(norm_count_files[[1]])
} else {
    stop(paste("There are", length(norm_count_files), "normalized count files in this directory. Please ensure there is one and try again."),
         call. = FALSE)
}
l2fc <-  data.table::fread(paste0(base_dir,"l2fc.csv")) ## ideally filename should be passed as arguments.
collapsed_l2fc <- data.table::fread(paste0(base_dir,"collapsed_l2fc.csv")) ## ideally filename should be passed as arguments.




## get negcon counts at day in QC table to pass to portal
cell_counts_negcon <- norm_counts %>% 
    dplyr::filter(!is.na(CCLE_name)) %>% 
    dplyr::filter(!(day %in% days_to_drop)) %>% 
    dplyr::filter(trt_type==control_type)

# validation: only one day should in QC table to pass to portal
if(length(cell_counts_negcon$day %>% unique())!=1) {
    print(paste("days_to drop arg was:", days_to_drop))
    print(paste("days left in df was:",cell_counts_negcon$day %>% unique()))
    stop("Error: should  have exactly one day in the QC table")
}
# validation: we need to have data in the dataframe
if(nrow(cell_counts_negcon) <1) {
    print(cell_counts_negcon)
    stop("Error: no negcon data at this day")
}
# validation: columns required in QC table should be present
if (!all(c("CCLE_name", "DepMap_ID", "cell_set", "day", "pert_plate") %in% colnames(cell_counts_negcon))){
    stop("All columns required in QC table not found, i.e:CCLE_name, DepMap_ID, cell_set, day, pert_plate")
}


## get median counts in negcon and apply count threshold to generate the QC table
med_cell_counts_qc <- cell_counts_negcon %>% 
    dplyr::group_by(CCLE_name, DepMap_ID, cell_set, day, pert_plate) %>% 
    dplyr::summarize(med_log_raw_counts=median(log2_n), 
                     med_log_norm_counts= median(log2_normalized_n),
                     med_raw_counts= median(n)) %>% 
    dplyr::mutate(pass_raw_count_qc=med_raw_counts>count_threshold) %>% 
    dplyr::ungroup() 
## threshold applied in raw counts where the pseudocount has not been added already, 
## in log2 transformed counts, the pseudocount has been added


## filter out the cell lines that failed QC from LFC output

collapsed_l2fc <- dplyr::left_join(collapsed_l2fc, 
                                   med_cell_counts_qc %>% 
                                       dplyr::select(CCLE_name, DepMap_ID, cell_set, day, pert_plate,
                                                     pass_raw_count_qc)) %>%
    dplyr::filter(pass_raw_count_qc==T) %>% 
    dplyr::select(-pass_raw_count_qc)

l2fc <- dplyr::left_join(l2fc, 
                         med_cell_counts_qc %>% 
                             dplyr::select(CCLE_name, DepMap_ID, cell_set, day, pert_plate, 
                                           pass_raw_count_qc)) %>%
    dplyr::filter(pass_raw_count_qc==T) %>% 
    dplyr::select(-pass_raw_count_qc)



## save the QC table output along with build name
out_path = paste0(out_dir, "/", build_name, "_EPS_QC_TABLE.csv")
print(paste0("Writing out EPS_QC table to ", out_path))
write.csv(med_cell_counts_qc, out_path,row.names= FALSE, quote= FALSE)


## overwrite the l2fc and collapsed_l2fc files with the filtered data
write.csv(l2fc, paste0(base_dir, "l2fc.csv"),row.names= FALSE, quote= FALSE)
write.csv(collapsed_l2fc, paste0(base_dir, "collapsed_l2fc.csv"), row.names= FALSE, quote= FALSE)










