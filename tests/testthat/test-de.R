metadata <- data.frame(
  batch = c("1", "2", "1", "2", "1", "2", "1", "2"),
  diagnosis = c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"),
  sex = c("M", "F", "M", "F", "M", "F", "M", "F"),
  RIN = c(5, 5, 5, 5, 5, 5, 5, 5),
  samples = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
  stringsAsFactors = FALSE
)

counts <- matrix(
  sample(0:100, size = 24),
  ncol = 8,
  dimnames = list(c("ENSG1", "ENSG2", "ENSG5"), c("S1", "S2", "S3", "S4",  "S5", "S6", "S7", "S8"))
)

metadata <- clean_covariates(metadata,
                             factors = c("samples", "sex", "batch", "diagnosis"),
                             continuous = c("RIN"),
                             sample_identifier = c("samples")
)

biomart_results <- tibble::column_to_rownames(
  tibble::tibble(ensembl_gene_id = rownames(counts), length = 400),
  var = "ensembl_gene_id"
)
de <- differential_expression(counts,
                              counts,
                              md = metadata,
                              primary_variable = c("sex", "diagnosis"),
                              biomart_results = biomart_results,
                              model_variables = c("batch"))

test_that("output is a list of 6", {
  expect_equal(length(de), 6)
})

test_that("voom output is complete", {
  expect_equal(names(de$voom_object), c("targets", "E", "weights"))
})

test_that("formula is complete", {
  expect_false(grepl("NULL", de$formula))
})
