biomart <- data.frame(
  row.names = c("ENSG1", "ENGS2", "ENSG3", "ENSG4"),
  chromosome_name = c("X", "X", "Y", "Y"),
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
  sample(0:100, size = 32),
  ncol = 8,
  dimnames = list(c("ENSG1", "ENSG2", "ENSG3", "ENSG4"),
                  c("S1", "S2", "S3", "S4",  "S5", "S6", "S7", "S8"))
))

plot <- plot_sexcheck_pca(metadata, counts, biomart, sex_var = "sex",
                          shape = "diagnosis", size = "RIN",
                          z = 2, split_condition = "sex")

test_that("output is plot", {
  expect_true("gg" %in% class(plot$plot))
})
