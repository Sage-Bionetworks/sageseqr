metadata <- data.frame(
  batch = rep( c("1", "2", "1", "2", "1", "2", "1", "2"), 2),
  diagnosis = rep( c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"), 2),
  sex = rep( c("M", "F", "M", "F", "M", "F", "M", "F"), 2),
  RIN = rep( c(5, 5, 5, 5, 5, 5, 5, 5), 2),
  row.names = c( paste0("S", c(1:16))),
  stringsAsFactors = TRUE
)

output <- build_formula(
  metadata,
  primary_variable = "diagnosis",
  model_variables = c("batch")
)

test_that("output data frame subset to primary and model variables.", {
  expect_equal(colnames(output$metadata), c("diagnosis", "batch"))
})

test_that("output is formula.",{
  expect_equal(class(output$formula), "formula")
})
