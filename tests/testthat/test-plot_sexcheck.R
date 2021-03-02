biomart <- data.frame(
  row.names = c("ENSG00000229807", "ENSG00000183878"),
  hgnc_symbol = c("XIST","UTY"),
  chromosome_name = c("X", "Y"),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  batch = c("1", "2", "1", "2", "1", "2", "1", "2"),
  diagnosis = c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"),
  sex = c("M", "F", "M", "F", "M", "F", "M", "F"),
  RIN = c(5, 5, 5, 5, 5, 5, 5, 5),
  row.names = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
  stringsAsFactors = FALSE
)

counts <- data.frame(matrix(
  sample(0:1000, size = 160),
  ncol = 80,
  dimnames = list(c("ENSG00000229807", "ENSG00000183878"),
                  c( paste0("S", c(1:80)))
)))

plot <- plot_sexcheck(metadata, counts, biomart_results = biomart, sex_var = "sex")

test_that("numeric values exist for each marker", {
  expect_true(sum(plot$sex_specific_counts$UTY) != 0)
  expect_true(sum(plot$sex_specific_counts$XIST) != 0)
})

test_that("there are no missing values, which would indicate a joining error", {
  expect_false(any(is.na(plot$sex_specific_counts$UTY)))
  expect_false(any(is.na(plot$sex_specific_counts$XIST)))
})

test_that("output is plot", {
  expect_true("gg" %in% class(plot$plot))
})
