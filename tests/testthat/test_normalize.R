#context("Test fastq read rawcounts")
library(magrittr)
library(prismSeqR)

test_that("Test normalize function", {
  filtered_counts = read.csv(test_path('assets', 'test_filtered_counts.csv'))
  CB_meta = read.csv(test_path('assets', 'CB_meta.csv'))
  
  normalized_counts = normalize(filtered_counts, CB_meta$Name )
  
  ground_truth = read.csv(test_path('assets', 'test_normalized_counts.csv'))
  
  expect_equal(('normalized_n' %in% colnames(normalized_counts)), TRUE)
  expect_equal(('log_normalized_n' %in% colnames(normalized_counts)), TRUE)
  expect_equal(nrow(ground_truth), nrow(normalized_counts))
}
)
