### If all data is correct format, test should pass ###
metadata <- data.frame(
  batch = c("1", "2", "1", "2"),
  diagnosis = c("dx", "dx", "ct", "ct"),
  row.names = c("S1", "S2", "S3", "S4"),
  stringsAsFactors = TRUE
)

counts <- data.frame(matrix(
  sample(0:100, size = 16),
  ncol = 4,
  dimnames = list(c("ENSG00000229807.12", "ENSG00000183878.12",
                    "ENSG00000239807.12", "ENSG00000259807.12"),
                  c("S1", "S2", "S3", "S4"))
))

test_that("passes if samples are same in both counts and metadata", {
  expect_null(check_mismatch(metadata, counts))
})

### Mismatch in samples between metadata and counts matrix ###
metadata <- data.frame(
  batch = c("1", "2", "1", "2"),
  diagnosis = c("dx", "dx", "ct", "ct"),
  row.names = c("S1", "S2", "S3", "S5"),
  stringsAsFactors = TRUE
)

counts <- data.frame(matrix(
  sample(0:100, size = 16),
  ncol = 4,
  dimnames = list(c("ENSG00000229807.12", "ENSG00000183878.12",
                    "ENSG00000239807.12", "ENSG00000259807.12"),
                  c("S1", "S2", "S3", "S4"))
))

test_that("errors if mismatched samples in metadata and counts matrix", {
  expect_error(check_mismatch(metadata, counts))
})
