suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))

#source("../src/load_libraries.R")
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readr)) #write_delim 
suppressPackageStartupMessages(library(stringr)) #str_detect
suppressPackageStartupMessages(library(dplyr)) #n(), %>%
suppressPackageStartupMessages(library(tidyr)) #pivot_wider
suppressPackageStartupMessages(library(reshape2))

## filter raw reads
## takes the raw readcount table and filters for expected indices and cell lines
## using the given metadata. QC metrics are written out to a text file
##
## takes:
##      raw_counts - an unfiltered counts table
##      sample_meta - the sample metadata for the particular experiment. Must follow the given set of guidelines for metadata
##      cell_line_meta - master metadata of cell lines
##      cell_set_meta - master metdata of cell sets and their contents
##      CB_meta - master metdata of control barcodes, their sequences, and their doses
filter_raw_reads = function(raw_counts, sample_meta, cell_line_meta, cell_set_meta, CB_meta, out=getwd()) {
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
  
  write(paste0("index_purity: ", round(index_purity*100, 2), "% \n",
               "cell_line_purity: ", round(cell_line_purity*100, 2), "%"),
        paste(out, "filtering_QC.txt", sep='/'))
  
  annotated_counts = cell_line_filtered %>% 
    merge(CB_meta, by.x="forward_read_cl_barcode", by.y="Sequence", all.x=T) %>% 
    dplyr::select_if(function(col) sum(is.na(col)) < length(col)) %>% 
    dplyr::select(-any_of(c("flowcell_name", "flowcell_lane", "index_1", "index_2", "members", 
                          "lysate_well", "lysate_plate", "pcr_well", "pcr_plate",
                          "forward_read_cl_barcode", "LUA"))) %>% 
    dplyr::relocate(project_code, CCLE_name, DepMap_ID, prism_cell_set, Name, log_dose, profile_id, sample_ID, trt_type, control_sample, control_barcodes,
                    bio_rep, tech_rep) %>% 
    dplyr::relocate(n, .after=last_col())
  
  return(annotated_counts)
}



## print_args
## writes configuration to file
##
## takes: 
##      args: args object from argparse
print_args <- function(args){
  config <- data.frame(args=names(args), values=unname(unlist(args)))
  config_path = paste(
    args$out, 
    "config.txt",
    sep="/"
  )
  print(paste("Saving config.txt file in :", config_path))
  write_delim(config, config_path, delim = ": ", col_names=F)
}

# create parser object
parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")

parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--raw_counts", default="raw_counts.csv", help = "path to file containing raw counts")
parser$add_argument("--out", default="", help = "Output path. Default is working directory")
parser$add_argument("-s", "--sample_meta", default="", help = "Sample metadata")
parser$add_argument("--cell_line_meta", default="../metadata/cell_line_meta.csv", help = "Cell Line metadata")
parser$add_argument("--cell_set_meta", default="../metadata/cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")
parser$add_argument("--id_cols", default="sample_ID,pcr_well,tech_rep", help = "Columns used to generate profile ids, comma-separated colnames from --sample_meta")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

#Save arguments
if (args$out == ""){
  args$out = args$wkdir
}
#print_args(args)

CB_meta = read.csv(args$CB_meta)
cell_line_meta = read.csv(args$cell_line_meta)
cell_set_meta = read.csv(args$cell_set_meta)
sample_meta = read.csv(args$sample_meta)
raw_counts = read.csv(args$raw_counts)

#split id_cols args
id_cols = unlist(strsplit(args$id_cols, ","))

if (!all(id_cols %in% colnames(sample_meta))){
  stop(paste("Colnames not found in sample_meta, check metadata or --id_cols argument:", args$id_cols))
}

#Generate Profile Ids based on sample_ID, PCR_well, and technical replicate 
#DONE: make profile_id columns parameters 
#sample_meta$profile_id = with(sample_meta, paste(sample_ID,pcr_well, tech_rep, sep = ":"))
sample_meta$profile_id = do.call(paste,c(sample_meta[id_cols], sep=':'))


print("creating filtered count file")
filtered_counts = filter_raw_reads(
  raw_counts,
  sample_meta,
  cell_line_meta,
  cell_set_meta,
  CB_meta,
  out=args$out
)

#Write Filtered Counts Table
filtrc_out_file = paste(
  args$out,
  'filtered_counts.csv',
  sep='/'
)

print(paste("writing filtered counts csv to: ", filtrc_out_file))
write.csv(filtered_counts, filtrc_out_file, row.names=F, quote=F)

## Make GCTX file
## For compatibility with other tools

# 
# counts_df = with(
#   filtered_counts[!is.na(filtered_counts$LUA),], #Remove rows where cell line LUA is not known
#   data.frame(
#     cid = profile_id,
#     rid=LUA,
#     value=n )
# )   #Extract matrix columns
# 
# #Pivot 
# counts_mat <- counts_df %>% 
#   pivot_wider(
#     names_from = cid,
#     id_cols = rid,
#     values_from = value,
#     values_fn=as.numeric
#   )
# 
# rids = counts_mat$rid
# counts_mat <- as.matrix(counts_mat[-1])
# rownames(counts_mat) <- rids
# 
# rownames(sample_meta) <- sample_meta$profile_id
# rownames(cell_line_meta) <- cell_line_meta$LUA
# 
# col_desc = sample_meta[sample_meta$profile_id %in% colnames(counts_mat),]
# row_desc = cell_line_meta[cell_line_meta$LUA %in% rids,]
# 
# #counts_gct <- new("GCT", mat=counts_mat, rid = rownames(counts_mat), cid = colnames(counts_mat), rdesc=row_desc, cdesc=col_desc)
# 
# filtrc_out_gct = paste(
#   args$out,
#   'Level2_counts',
#   sep='/'
# )
# #print(paste0("Writing GCTx file to ", filtrc_out_gct))
# #write_gctx(counts_gct, filtrc_out_gct)

