library(argparse)
library(prismSeqR)
library(magrittr)

parser <- ArgumentParser()
# specify our desired options 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-c", "--filtered_counts", default="filtered_counts.csv",
                    help="path to file containing filtered counts")
parser$add_argument("--CB_meta", default="../metadata/CB_meta.csv", help = "Control Barcode metadata")
parser$add_argument("-o", "--out", default=getwd(), help = "Output path. Default is working directory")
parser$add_argument("--pseudocount", default=20, help = "pseudo count for normalization")

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

filtered_counts = read.csv(args$filtered_counts)
CB_meta = read.csv(args$CB_meta)
pseudocount_arg= args$pseudocount
pseudocount = as.numeric(pseudocount_arg)
print("creating normalized count file")
normalized_counts = filtered_counts %>% 
  normalize(CB_meta$Name,pseudocount) 

normcounts_out_file = paste(
  args$out,
  "normalized_counts.csv",
  sep='/'
)

if (args$verbose){
  print(paste("Writing normalized count file", normcounts_out_file))
}

write.csv(normalized_counts, normcounts_out_file, row.names=F, quote=F)


