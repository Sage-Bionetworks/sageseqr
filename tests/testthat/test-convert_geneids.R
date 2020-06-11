test_that("No warning is returned", {
  df <- tibble::tribble(~genes, ~sample1,
                "ENSG400", 4,
                "ENSG401", 6,
                "ENSG402", 10) %>%
    tibble::column_to_rownames(., var = "genes")
  testthat::expect_silent(convert_geneids(df))
})
