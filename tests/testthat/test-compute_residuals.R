metadata <- data.frame(
  batch = c("1", "2", "1", "2", "1", "2", "1", "2"),
  diagnosis = c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"),
  sex = c("M", "F", "F", "F", "M", "F", "F", "F"),
  RIN = c(5, 5, 4, 5, 4, 4, 5, 5),
  site = c("one", "one", "one", "two", "two", "three", "three", "three"),
  row.names = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
  stringsAsFactors = TRUE
)

counts <- data.frame(matrix(
  sample(0:100, size = 32),
  ncol = 8,
  dimnames = list(c("ENSG00000229807.12", "ENSG00000183878.12",
                    "ENSG00000239807.12", "ENSG00000259807.12"),
                  c("S1", "S2", "S3", "S4",  "S5", "S6", "S7", "S8"))
))

matrix <- compute_residuals(
  metadata,
  counts,
  counts,
  primary_variable = c("diagnosis", "sex"),
)

test_that("gene features stored in column", {
  expect_true("feature" %in% colnames(matrix$output))
})

test_that("predictor variable preserved as signal", {
  expect_true(all(grepl("diagnosis_sex", matrix$signal)))
})

test_that("output list meets expectations", {
  expect_identical(names(matrix), c("output", "signal", "adjusted_fit", "formula"))
})
