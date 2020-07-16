metadata <- tibble::tribble(
  ~samples, ~sex, ~batch, ~age, ~RIN, ~diagnosis,
  "S567", "F", "1", "68", "7.5", "dx",
  "S453", "F", "2", "60", "7", "dx",
  "S231", "M", "2", "54", "7", "dx",
  "S444", "M", "2", "55", "7", "dx"
)

metadata <- tibble::column_to_rownames(metadata, var = "samples")

counts <- tibble::tribble(
  ~geneId, ~S567, ~S453, ~S444, ~S231,
  "ENSG1", 0, 0, 500, 0,
  "ENSG2", 23, 25, 29, 30,
  "ENSG3", 84, 82, 0, 0,
  "ENSG4", 0, 0, 0, 0,
  "ENSG5", 0, 0, 0, 70
)

counts <- tibble::column_to_rownames(counts, var = "geneId")

test_that("zero count genes are removed", {
  filtered <- filter_genes(metadata, counts, conditions = "sex",
                           cpm_threshold = 1, conditions_threshold = 0.5)
  expect_equal(c("ENSG1", "ENSG2", "ENSG3", "ENSG5"), rownames(filtered))
})
