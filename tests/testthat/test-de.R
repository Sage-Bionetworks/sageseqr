library(variancePartition)
data(varPartDEdata)

countMatrix <- countMatrix[1:100,]

# Add continuous variable
metadata$AlignmentRate <- 0.9

rownames(metadata) <- NULL

metadata <- clean_covariates(metadata,
                             factors = c("Individual", "Disease", "Experiment", "Sex"),
                             continuous = c("AlignmentRate"),
                             sample_identifier = c("Experiment"))

biomart_results <- tibble::column_to_rownames(
  tibble::tibble(ensembl_gene_id = rownames(countMatrix), length = 400),
  var = "ensembl_gene_id"
  )

test_that("output is a list of 4", {
  de <- differential_expression(countMatrix,
                                countMatrix,
                                md = metadata,
                                model_variables = c("Individual"),
                                primary_variable = c("Sex", "Disease"),
                                biomart_results = biomart_results)
  expect_equal(length(de), 4)
})
