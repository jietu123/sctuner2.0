
packages <- commandArgs(trailingOnly = TRUE)
for (pkg in packages) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  version <- ""
  load_success <- FALSE
  error_message <- ""
  if (installed) {
    version <- as.character(utils::packageVersion(pkg))
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      load_success <- TRUE
    }, error = function(e) {
      error_message <<- conditionMessage(e)
    })
  } else {
    error_message <- "package not installed"
  }
  cat(pkg, installed, version, load_success, gsub("\t|\n", " ", error_message), sep="\t")
  cat("\n")
}
