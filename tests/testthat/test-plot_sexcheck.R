biomart <- data.frame(
  row.names = c("ENSG00000229807", "ENSG00000183878"),
  hgnc_symbol = c("XIST","UTY"),
  chromosome_name = c("X", "Y"),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  batch = rep( c("1", "2", "1", "2", "1", "2", "1", "2"), 10),
  diagnosis = rep( c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"), 10),
  sex = rep( c("M", "F", "M", "F", "M", "F", "M", "F"), 10),
  RIN = rep( c(5, 5, 5, 5, 5, 5, 5, 5), 10),
  row.names = c( paste0("S", c(1:80))),
  stringsAsFactors = FALSE
)

counts <- data.frame(matrix(
  sample(0:1000, size = 160),
  ncol = 80,
  dimnames = list(c("ENSG00000229807", "ENSG00000183878"),
                  c( paste0("S", c(1:80)))
)))

plot <- plot_sexcheck(clean_metadata=metadata, count_df=counts, biomart_results = biomart, sex_var = "sex")

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
