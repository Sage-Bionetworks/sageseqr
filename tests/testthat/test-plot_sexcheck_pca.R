biomart <- data.frame(
  row.names = c("ENSG00000229807", "ENSG00000183878",
                "ENSG00000183811", "ENSG00000229822",
                "ENSG00000183833", "ENSG00000229844",
                "ENSG00000XXXX", "ENSG00000ZZZZ"),
  hgnc_symbol = c("XIST","UTY",
                  "Foo1", "Foo2",
                  "Foo3", "Foo4",
                  "Foo5", "Foo6"
                  ),
  chromosome_name = c("X", "Y",
                      "X", "Y",
                      "X", "Y",
                      "X", "Y"),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  batch = rep( c("1", "2", "1", "2", "1", "2", "1", "2"), 10),
  diagnosis = c(rep("dx", 20), rep("ct", 20),rep("dx", 20), rep("ct", 20)),
  sex = c(rep("M", 40), rep("F", 40)),
  RIN = rep( c(5, 5, 5, 5, 5, 5, 5, 5), 10),
  row.names = c( paste0("S", c(1:80))),
  stringsAsFactors = FALSE
)

ENSG00000183878 <- c(rchisq(40,24), rchisq(40,2))
ENSG00000229807 <- c(rchisq(40,2), rchisq(40,24))
ENSG00000183811 < -c(rchisq(53,5), rchisq(27,4))
ENSG00000229822 <- c(rchisq(53,3), rchisq(27,7))
ENSG00000183833 <- c(rchisq(17,16), rchisq(63,22))
ENSG00000229844 <- c(rchisq(17,98), rchisq(63,97))
ENSG00000XXXX <- c(rchisq(20,80), rchisq(20,5), rchisq(20,80), rchisq(20,5))
ENSG00000ZZZZ <- c(rchisq(20,80), rchisq(20,5), rchisq(20,80), rchisq(20,5))
counts <- as.data.frame(rbind(ENSG00000229807,ENSG00000183878,
                               ENSG00000183811, ENSG00000229822,
                               ENSG00000183833, ENSG00000229844,
                               ENSG00000XXXX, ENSG00000ZZZZ))
colnames(counts) <- paste0( 'S', 1:80)

plot <- plot_sexcheck_pca(metadata, counts, biomart, sex_var = "sex")

test_that("output is plot", {
  expect_true("gg" %in% class(plot$plot))
})
