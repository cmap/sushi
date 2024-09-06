#context("Test filter counts")
source("filter_counts.R")
library(magrittr)

test_that("Test filter counts function", {
  expect_error(validate_cell_line_meta("Index1"), "Cell line meta validation not yet implemented")
}
)
test_that("Test LUA duplicates in cell line meta with duplicates", {
  cell_line_meta_with_duplicates <- data.frame(LUA = c("LUA-1000", "LUA-1001", "LUA-1002", "LUA-1001"))
  expected_result <- data.frame(LUA = c("LUA-1000", "LUA-1001", "LUA-1002"))
  luas_removed = remove_duplicate_luas(cell_line_meta_with_duplicates)
  expect_equal(luas_removed$LUA %>% sort(), expected_result$LUA %>% sort()) 
  expect_equal(length(luas_removed$LUA), length(expected_result$LUA)) 
  })

# Classes are different between LUAs removed and expected result - LUA column needs to be subsetted for comparison
test_that("Test LUA duplicates in cell line meta without duplicates", {
  cell_line_meta <- data.frame(LUA = c("LUA-1000", "LUA-1001", "LUA-1002"))
  expected_result <- data.frame(LUA = c("LUA-1000", "LUA-1001", "LUA-1002"))
  luas_removed = remove_duplicate_luas(cell_line_meta)
  expect_equal(luas_removed$LUA, expected_result$LUA) 
})


test_that("Test merge assay pool meta function with pool_id", {
  test_assay_pool_meta <- data.frame(depmap_id = "ACH-000956", ccle_name = "22RV1_PROSTATE",
                                     davepool_id = "EXT.PR500.CS01.2.B", pool_id = "P111")
  
  test_filtered_counts <- data.frame(CCLE_name = "22RV1_PROSTATE", cell_set = "EXT.PR500.CS01.2.B", 
                                     DepMap_ID = "ACH-000956", random_column = "X")
  
  davepool_ids = c("EXT.PR500.CS01.2.B")
  expected_result = test_filtered_counts
  expected_result$pool_id = "P111"
  expect_equal(merge_pool_info(test_filtered_counts, test_assay_pool_meta, davepool_ids), expected_result)
})

test_that("Test merge assay pool meta function with no pool_id", {
  test_assay_pool_meta <- data.frame(depmap_id = "ACH-000956", ccle_name = "22RV1_PROSTATE",
                                     davepool_id = "EXT.PR500.CS01.2.B")
  
  test_filtered_counts <- data.frame(CCLE_name = "22RV1_PROSTATE", cell_set = "EXT.PR500.CS01.2.B", 
                                     DepMap_ID = "ACH-000956", random_column = "X")
  
  davepool_ids = c("EXT.PR500.CS01.2.B")
  expect_error(merge_pool_info(test_filtered_counts, test_assay_pool_meta, davepool_ids), "Column `pool_id` doesn't exist.")
})

