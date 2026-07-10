
suppressPackageStartupMessages(library(jsonlite))
payload <- list(
  R_version = R.version.string,
  R_home = R.home(),
  libPaths = .libPaths(),
  locale = Sys.getlocale(),
  encoding = getOption("encoding"),
  platform = R.version$platform
)
cat(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE))
