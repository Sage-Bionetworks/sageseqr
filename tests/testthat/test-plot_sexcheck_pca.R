biomart <- data.frame(
  row.names = c("ENSG00000229807", "ENSG00000183878"),
  hgnc_symbol = c("XIST","UTY"),
  chromosome_name = c("X", "Y"),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  batch = rep( c("1", "2", "1", "2", "1", "2", "1", "2"), 10),
  diagnosis = rep( c("dx", "dx", "ct", "ct", "dx", "dx", "ct", "ct"), 10),
  sex = rep( c("M", "F", "M", "F", "M", "F", "M", "F"), 10),
  RIN = rep( c(5, 5, 5, 5, 5, 5, 5, 5), 10),
  row.names = c( paste0("S", c(1:80))),
  stringsAsFactors = FALSE
)

counts <- data.frame(matrix(
  c(sample(0:1000, size = 160), sample(750:1200, size = 160), sample(1750:2200, size = 160), sample(550:625, size = 160, replace = TRUE)),
  ncol = 80,
  dimnames = list(c("ENSG00000229807", "ENSG00000183878", "ENSG00000XXXX", "ENSG00000YYYY", "ENSG00000ZZZZ", "ENSG00000AAAA", "ENSG00000BBBB", "ENSG00000CCCC"),
                  c( paste0("S", c(1:80)))
)))

plot <- plot_sexcheck_pca(clean_metadata = metadata, filtered_counts = counts, biomart_results = biomart, sex_var = "sex")

test_that("output is plot", {
  expect_true("NULL" %in% class(plot$plot))
})
