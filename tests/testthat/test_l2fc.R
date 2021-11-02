#context("Test fastq read rawcounts")
library(magrittr)
library(prismSeqR)

test_that("Test l2fc function", {
  
  normalized_counts = read.csv(test_path('assets', 'test_normalized_counts.csv'))
  l2fc_table = compute_l2fc(normalized_counts, control_type = "negcon")
  
  ground_truth = read.csv(test_path('assets', 'test_l2fc.csv'))
  
  expect_equal(('l2fc' %in% colnames(l2fc_table)), TRUE)
  expect_equal(nrow(ground_truth), nrow(l2fc_table))
  
}
)
