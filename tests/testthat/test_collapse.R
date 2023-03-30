#context("Test collapse function")
library(magrittr)
library(prismSeqR)

test_that("Test collapse function", {
  
  l2fc_table = read.csv(test_path('assets', 'test_l2fc.csv'))
  l2fc_collapsed = collapse_counts(l2fc_table %>% dplyr::mutate(mean_normalized_n = 0))
  
  ground_truth = read.csv(test_path('assets', 'test_l2fc_collapsed.csv'))
  
  expect_equal(('median_l2fc' %in% colnames(l2fc_collapsed)), TRUE)
  expect_equal(nrow(ground_truth), nrow(l2fc_collapsed))
}
)