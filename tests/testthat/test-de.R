metadata <- tibble::tribble(
  ~samples, ~sex, ~batch, ~age, ~RIN, ~diagnosis,
  "S567", "F", "1", "68", "7.5", "dx",
  "S453", "F", "2", "60", "7", "dx",
  "S231", "M", "2", "54", "7", "dx",
  "S444", "M", "2", "55", "7", "dx",
  "S555", "F", "1", "68", "7.5", "dx",
  "S304", "F", "2", "60", "7", "ct",
  "S232", "M", "2", "54", "7.5", "ct",
  "S447", "M", "2", "55", "7.5", "ct"
)

counts <- tibble::tribble(
  ~geneId, ~S567, ~S453,~S231, ~S444, ~S555, ~S304, ~S232, ~S447,
  "ENSG1", 0, 0, 500, 0, 10, 43, 432, 45,
  "ENSG2", 23, 25, 29, 30, 56, 89, 73, 70,
  "ENSG3", 84, 82, 0, 0, 98, 80, 200, 201,
  "ENSG4", 0, 204, 350, 790, 353, 456, 555, 890,
  "ENSG5", 0, 0, 0, 70, 0, 0, 0, 70
)

counts <- tibble::column_to_rownames(counts, var = "geneId")

counts <- as.matrix(counts)

metadata <- clean_covariates(metadata,
                             factors = c("samples", "sex", "batch", "diagnosis"),
                             continuous = c("RIN", "age"),
                             sample_identifier = c("samples")
)

biomart_results <- tibble::column_to_rownames(
  tibble::tibble(ensembl_gene_id = rownames(counts), length = 400),
  var = "ensembl_gene_id"
)
de <- differential_expression(counts,
                              counts,
                              md = metadata,
                              model_variables = c("batch"),
                              primary_variable = c("sex", "diagnosis"),
                              biomart_results = biomart_results)

test_that("output is a list of 5", {
  expect_equal(length(de), 5)
})

