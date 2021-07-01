metadata <- data.frame(
  batch = c("1", "2", "1", "2", "1", "2", "1", "2"),
  diagnosis = c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"),
  sex = c("M", "F", "M", "F", "M", "F", "M", "F"),
  RIN = c(5, 5, 5, 5, 5, 5, 5, 5),
  sampleID = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
  stringsAsFactors = FALSE
)

counts <- data.frame(matrix(
  sample(seq(0, 100, by = 0.01), size = 16),
  ncol = 8,
  dimnames = list(c("ENSG00000229807.12", "ENSG00000183878.12"),
                  c("S1", "S2", "S3", "S4",  "S5", "S6", "S7", "S8"))
))

counts_df <- tibble::rownames_to_column(counts, "feature")

test_that("function runs and computation is correct", {
  dat <- compute_mean_sd(metadata, "sampleID", counts_df, "feature")
  #drop gene feature to test row means and row sds
  test <- counts_df[,-1]
  testthat::expect_equal(colnames(dat), c("feature", "n", "mean", "sd"))
  testthat::expect_equal(sd(test[1,], na.rm = TRUE), dat$sd[1])
  testthat::expect_equal(rowMeans(test, na.rm = TRUE), dat$mean)
  testthat::expect_equal(dim(metadata)[1], 8)
})
