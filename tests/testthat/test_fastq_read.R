#context("Test fastq read rawcounts")
library(magrittr)
library(prismSeqR)

test_that("fastq read of small datasets finds 12 reads", {
  fastq_test_files = test_path('assets', 'fastq')  
  
  read_directory_contents <- fastq_test_files %>% 
    purrr::map(list.files, full.names = T) %>%
    purrr::reduce(union)
  
  barcode_read_files <- read_directory_contents %>%
    purrr::keep(stringr::str_detect, fixed('_R1_')) %>%
    sort()
  
  index_1_files <- read_directory_contents %>%
    purrr::keep(stringr::str_detect, fixed('_I1_')) %>%
    sort()
  
  index_2_files <- read_directory_contents %>%
    purrr::keep(stringr::str_detect, fixed('_I2_')) %>%
    sort()
  
  
  raw_counts <- write_df_from_fastq(forward_read_fastq_files = barcode_read_files,
                                    index_1_file = index_1_files,
                                    index_2_file = index_2_files)
  
  expect_equal(sum(raw_counts$n), 12)
  expect_equal(length(raw_counts), 4)
  expect_equal(dim(raw_counts)[1], 5)
}
)