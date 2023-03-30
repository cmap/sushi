#context("Test filter counts")
library(magrittr)
library(prismSeqR)

test_that("Test filter counts function", {
  test_dir = test_path('assets')
  
  rc = read.csv(paste(test_dir, 'test_raw_counts.csv',sep='/'))
  
  sm = read.csv(paste(test_dir, 'test_meta.csv',sep='/'))
  
  cell_meta = read.csv(paste(test_dir, 'test_cell_line_meta.csv',sep='/'))
  cell_set_meta = read.csv(paste(test_dir, 'test_cell_set_meta.csv',sep='/'))
  CB_meta = read.csv(paste(test_dir, 'CB_meta.csv',sep='/'))
  
  filtered_counts = filter_raw_reads(rc, sm, cell_meta, cell_set_meta, CB_meta, id_cols = c('treatment', 'dose','dose_unit','day') )
  
  expect_equal(nrow(filtered_counts$filtered_counts), 69)
  expect_equal(length(filtered_counts$filtered_counts), 17)
}
)