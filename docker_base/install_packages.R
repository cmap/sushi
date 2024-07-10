options(repos = c(
  "https://iwww.broadinstitute.org/~datasci/R-packages",
  "https://cran.rstudio.com/"))

options(repos=structure(c(CRAN="http://cran.r-project.org")))

install_and_log <- function(pkg) {
  tryCatch(
    {
      install.packages(pkg)
      message(sprintf("Successfully installed %s", pkg))
    },
    error = function(e) {
      message(sprintf("Failed to install %s: %s", pkg, e$message))
    }
  )
}

bioc_install_and_log <- function(pkg) {
  tryCatch(
    {
      BiocManager::install(pkg, ask = FALSE)
      message(sprintf("Successfully installed %s", pkg))
    },
    error = function(e) {
      message(sprintf("Failed to install %s: %s", pkg, e$message))
    }
  )
}

# Install dependencies for tidyverse
install_and_log("googledrive")
install_and_log("googlesheets4")
install_and_log("haven")
install_and_log("httr")
install_and_log("ragg")
install_and_log("rvest")
install_and_log("xml2")

# Install tidyverse
install_and_log("tidyverse")

install_and_log("hdf5r")
install_and_log("reshape2")
install_and_log("readr")
install_and_log("magrittr")
install_and_log("dr4pl")
install_and_log("data.table")
install_and_log("scam")
install_and_log("argparse")
install_and_log("splitstackshape")
install_and_log("BiocManager")
BiocManager::install()

bioc_install_and_log("sva")
bioc_install_and_log("limma")
bioc_install_and_log("apeglm")
bioc_install_and_log("ShortRead")

install_and_log("PRROC")
install_and_log("gmodels")
install_and_log("R.utils")

# Uncomment these lines if you need the devtools packages
#install_and_log("devtools")
#tryCatch(
#  devtools::install_github("broadinstitute/cdsr_models", dependencies=TRUE, force=TRUE),
#  error = function(e) message(sprintf("Failed to install cdsr_models: %s", e$message))
#)
#tryCatch(
#  devtools::install_github("broadinstitute/cdsr_biomarker"),
#  error = function(e) message(sprintf("Failed to install cdsr_biomarker: %s", e$message))
#)
#tryCatch(
#  devtools::install_github('https://github.com/cmap/sushi'),
#  error = function(e) message(sprintf("Failed to install sushi: %s", e$message))
#)

