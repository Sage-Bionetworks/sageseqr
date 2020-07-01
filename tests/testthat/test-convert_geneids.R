test_that("No warning is returned", {
  df <- tibble::tribble(~genes, ~sample1,
                "ENSG400", 4,
                "ENSG401", 6,
                "ENSG402", 10) %>%
    tibble::column_to_rownames(., var = "genes")
  testthat::expect_silent(convert_geneids(df))
})

test_that("Transcript Ids removed", {
  expected_output <- c("ENSG400", "ENSG401", "ENSG402")
  df <- tibble::tribble(~genes, ~sample1,
                        "ENSG400.5", 4,
                        "ENSG401.7", 6,
                        "ENSG402.8", 10) %>%
    tibble::column_to_rownames(., var = "genes")
  output <- convert_geneids(df)
  testthat::expect_equal(output, expected_output)
})
