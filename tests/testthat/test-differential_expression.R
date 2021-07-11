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

biomart_results <- data.frame(
  gene_length = c(200,200),
  row.names = c("ENSG00000229807.12", "ENSG00000183878.12")
)

test_that("differential_expression executes completely", {
  expect_length(
    differential_expression(
      counts,
      counts,
      metadata,
      "diagnosis",
      biomart_results,
      0.05,
      1.5
      ),
    6
  )
})
