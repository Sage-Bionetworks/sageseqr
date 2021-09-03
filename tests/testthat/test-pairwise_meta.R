dat <- tibble::tribble(~tissue, ~diagnosis, ~sex, ~feature, ~n, ~mean, ~sd,
               "cbe", "dx", "f", "ENSG00000225630", 47, -0.24, 1.256,
               "pf", "dx", "m", "ENSG00000237973", 54, -0.23, 0.1,
               "cbe", "dx", "m", "ENSG00000237973", 47, -0.4, 0.213,
               "cbe", "dx", "f", "ENSG00000237973", 47, -0.277, 0.56,
               "cbe", "other", "m", "ENSG00000237973", 72,1.003, 0.5,
               "cbe", "ct", "f", "ENSG00000225630", 88, -2.0, 2.2,
               "cbe", "ct", "f", "ENSG00000237973", 88, -1.6, 0.44,
               "cbe", "other", "f", "ENSG00000225630", 72, 1.34, 1.01,
               "pf", "other", "m", "ENSG00000237973", 66,0.98, 0.7,
               "cbe", "ct", "m", "ENSG00000237973", 88,-1, 0.24,
               "pf", "dx", "f", "ENSG00000237973", 54, -0.16, 0.5,
               "pf", "dx", "f", "ENSG00000225630", 54, -0.3, 1.1,
               "pf", "ct", "m", "ENSG00000237973", 57,-0.87, 0.2,
               "pf", "other", "f", "ENSG00000225630", 66, 1.03, 0.98,
               "pf", "other", "f", "ENSG00000237973", 66, 1.6, 0.8,
               "pf", "ct", "f", "ENSG00000225630", 57, -2.06, 2.0,
               "cbe", "other", "f", "ENSG00000237973", 72, 1.78, 0.7,
               "pf", "ct", "f", "ENSG00000237973", 57, -1.67, 0.3
)
# pairwise_meta() requires a constrained format where summary statistics are
# grouped by a condition of interest
input <- dplyr::group_by(dat, diagnosis, sex, feature) %>%
  tidyr::nest() %>%
  dplyr::arrange(feature, .by_group = TRUE)

index <- input[,c("sex", "feature")] %>%
  dplyr::distinct()

meta_stats <- pairwise_meta(c("dx", "ct"), input, "diagnosis", index)

test_that("expect dimensions of output to be 3", {
  expect_equal(3, dim(meta_stats)[1])
})

test_that("columns of output match expectations", {
  expect_equal(colnames(meta_stats)[1:2], c("sex", "feature"))
})

test_that("fails when more than two comparison_values are passed",{
  expect_error(
    pairwise_meta(c("dx", "ct", "ct"), input, diagnosis, index)
  )
})
