suppressPackageStartupMessages(library(argparse))
#suppressMessages(library(cmapR))
suppressPackageStartupMessages(library(dplyr)) #n()
suppressPackageStartupMessages(library(scam))
suppressPackageStartupMessages(library(magrittr))
#suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(grDevices))

## filteredCounts_QC
## checks that technical and biological replicates of each sample are all well correlated with each other
## checks that control barcodes follow expected linear trend in all samples
## saves images of sample correlation and control barcode trends
## 
## takes:
##      filtered_counts - dataframe of filtered counts 
##      out - the path to the folder where QC images should be saved
filteredCounts_QC = function(filtered_counts, cell_set_meta, out) {
  num_profiles = filtered_counts$profile_id %>% unique() %>% length()
  
  # counts in each sample, colored by cell line or control barcode
  total_counts = filtered_counts %>% ungroup() %>% 
    mutate(type = ifelse(is.na(Name), "cell line", "control barcode")) %>% 
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
  
  # fraction of cexpected cell lines in each sample
  num_cell_lines = filtered_counts %>% ungroup() %>% 
    filter(is.na(Name)) %>% 
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
  correlation_matrix = filtered_counts %>% ungroup() %>% 
    filter(is.na(Name)) %>% 
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
  cbt = filtered_counts %>% ungroup() %>% 
    filter(is.na(CCLE_name)) %>% 
    ggplot(aes(x=log_dose, y=log10(n))) +
    geom_point() +
    geom_smooth(method = 'glm') +
    stat_regline_equation(aes(label =  ..eq.label..), label.y.npc = "top") +
    stat_regline_equation(aes(label =  ..adj.rr.label..), label.y.npc = "top", position = position_nudge(y=-0.5)) +
    facet_wrap(~profile_id) +
    labs(x="log(dose)")
  
  pdf(file=paste(out, "control_barcode_trend.pdf", sep="/"),
      width=sqrt(num_profiles)*2, height=sqrt(num_profiles)*2)
  print(cbt)
  dev.off()
  
  # cell line counts vs. control barcode counts
  # assumes that last piece of profile_id is tech_rep
  aclt = filtered_counts %>% ungroup() %>% 
    mutate(type = ifelse(is.na(Name), "cell line", "control barcode"),
           log_n = log10(n),
           sample_id = gsub('.{2}$', '', profile_id)) %>% 
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


parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("--cell_set_meta", default="../metadata/cell_set_meta.csv", help = "Cell set metadata")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

filtered_counts = read.csv(args$filtered_counts)
cell_set_meta = read.csv(args$cell_set_meta)

print("generating filtered counts QC images")
filteredCounts_QC(filtered_counts, cell_set_meta, args$out)
