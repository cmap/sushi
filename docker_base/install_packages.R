install.packages("argparse")
install.packages("magrittr")
install.packages("tidyverse")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cmapR")