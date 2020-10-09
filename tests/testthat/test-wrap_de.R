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

conditions <- list(
  interaction = c("diagnosis", "sex"),
  sex = c("sex"),
  diagnosis = c("diagnosis")
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

all_de <- wrap_de(conditions,
                  counts,
                  counts,
                  metadata,
                  model_variables = c("batch"),
                  biomart_results)

test_that("conditions match differential expression output", {
  expect_equal(names(conditions), names(all_de))
})

test_that("all conditions execute", {
  expect_equal(length(all_de), length(conditions))
})

test_that("differential expression data frame exists", {
  expect_true(is.data.frame(all_de$interaction$differential_expression))
})
