test_that("sapply catches rownames", {
  dat <- tibble::tribble(~gene_ids, ~samp1, ~samp2,
                         "N_", 456, 6789,
                         "N_something", 43234, 54643,
                         "ENSG1111", 0, 0)
  dat <- tibble::column_to_rownames(dat, var = "gene_ids")

  dat <- parse_counts(dat)

  expect <- tibble::tribble(~gene_ids, ~samp1, ~samp2,
                            "ENSG1111", 0, 0) %>%
    tibble::column_to_rownames(., var = "gene_ids")
  expect_equal(expect, dat)
})
