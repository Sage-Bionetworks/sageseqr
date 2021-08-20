library(targets)
library(config)
  list(
    tar_target(
      compute_model,
      c("a list", "of", "vars"),
    ),
    tar_target(
      model,
      if(is.null(get("de log2F"))) {compute_model} else {get("de log2F")}
    ),
    tar_target(
      store_output,
      write.csv(model, "test.csv")
    )
  )
