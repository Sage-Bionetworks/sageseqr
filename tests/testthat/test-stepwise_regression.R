metadata <- data.frame(
    batch = c("1", "2", "1", "2", "1", "2", "1", "2"),
    diagnosis = c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"),
    sex = c("M", "F", "M", "F", "M", "F", "M", "F"),
    RIN = c(5, 5, 5, 5, 5, 5, 5, 5),
    row.names = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"),
    stringsAsFactors = FALSE
  )
  
counts <- data.frame(matrix(
    sample(0:100, size = 16),
    ncol = 8,
    dimnames = list(c("ENSG00000229807.12", "ENSG00000183878.12"),
                    c("S1", "S2", "S3", "S4",  "S5", "S6", "S7", "S8"))
  ))

test_that("Skip Function in Stepwise Model Regression is implemented", {
  expect_identical(stepwise_regression( metadata, 
                              model_variables = c('batch','RIN'), 
                              primary_variable="diagnosis",
                              counts,
                              skip = TRUE),
         "Skipping Stepwise Model Generation")
})
