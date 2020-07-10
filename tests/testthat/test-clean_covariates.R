metadata <- tibble::tribble(
  ~samples, ~sex, ~batch, ~age, ~RIN,
  "S567", "F", "1", "68", "7.5",
  "S453", "F", "1", "60", "7",
  "S231", "M", "2", "54", "7"
)

factor_variables <- c("sex", "batch")
continuous_variables <- c("age", "RIN")

test_that("factors are factors", {
  md <- clean_covariates(metadata,
                         factors = factor_variables,
                         continuous = continuous_variables,
                         sample_identifier = "samples")
  expect_true(all(sapply(md[,factor_variables], is.factor)))
})

test_that("numeric variables are numeric", {
  md <- clean_covariates(metadata,
                         factors = factor_variables,
                         continuous = continuous_variables,
                         sample_identifier = "samples")
  expect_true(all(sapply(md[,continuous_variables], is.numeric)))
})

df_metadata <- as.data.frame(metadata)

test_that("function accepts data frames", {
  md <- clean_covariates(df_metadata,
                   factors = factor_variables,
                   continuous = continuous_variables,
                   sample_identifier = "samples")
  expect_true(all(sapply(md[,factor_variables], is.factor)))
  expect_true(all(sapply(md[,continuous_variables], is.numeric)))
})

test_that("user cannot omit factor and continuous vector arguments", {
  expect_error(clean_covariates(df_metadata))
})

test_that("factor and continous variables are discrete.",{
  expect_error(clean_covariates(df_metadata,
                                factors = c("sex"),
                                continuous = c("sex"),
                                sample_identifier = "samples"))
})

test_that("variables are in md", {
  expect_error(clean_covariates(metadata,
                                factors = c("sex"),
                                continuous = c("diagnosis"),
                                sample_identifier = "samples"))
})

test_that("rownames are equivalent to the sample identifiers provided", {
  df_metadata <- clean_covariates(metadata,
                                  factors = factor_variables,
                                  continuous = continuous_variables,
                                  sample_identifier = "samples")
  expect_equal(rownames(df_metadata), metadata$samples)
})
