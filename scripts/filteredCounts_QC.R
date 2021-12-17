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
filteredCounts_QC = function(filtered_counts, out) {
  correlation_matrix = filtered_counts %>% 
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
      width=10, height=10)
  print(cp)
  dev.off()
  
  cbt = filtered_counts %>% 
    filter(is.na(CCLE_name)) %>% 
    ggplot(aes(x=log_dose, y=log10(n))) +
    geom_point() +
    geom_smooth(method = 'glm') +
    stat_regline_equation(aes(label =  ..eq.label..), label.y.npc = "top") +
    stat_regline_equation(aes(label =  ..adj.rr.label..), label.y.npc = "top", position = position_nudge(y=-0.5)) +
    facet_wrap(~profile_id) +
    labs(x="log(dose)")
  
  pdf(file=paste(out, "control_barcode_trend.pdf", sep="/"),
      width=10, height=10)
  print(cbt)
  dev.off()
}


parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("--wkdir", default=getwd(), help="Working directory")
parser$add_argument("-c", "--normalized_counts", default="normalized_counts.csv",
                    help="path to file containing normalized counts")
parser$add_argument("-o","--out", default="", help = "Output path. Default is working directory")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

if (args$out == ""){
  args$out = args$wkdir
}

normalized_counts = read.csv(args$normalized_counts)

print("checking replicate correlation")
check_replicate_cor(normalized_counts, args$out)
