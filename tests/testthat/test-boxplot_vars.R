test_that("output is plot", {
  metadata <- tibble::tribble(~diagnosis, ~ageofdeath, ~pmi,
                              "CT", 60, 4,
                              "CT", 61, 5,
                              "ZZ", 50, 4)
  output <- boxplot_vars(metadata,
                         include_vars = c("pmi", "ageofdeath"),
                         x_var = "diagnosis")
  testthat::expect_equal(class(output), c("gg", "ggplot"))
})
