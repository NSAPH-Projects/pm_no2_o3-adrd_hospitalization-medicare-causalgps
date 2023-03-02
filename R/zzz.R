.onLoad <- function(libname, pkgname) {

  pb_internal <- "study_data/public/internal"
  if (!dir.exists(pb_internal)) {
    dir.create(pb_internal, recursive = TRUE)
  }

  pr_internal <- "study_data/private/internal"
  if (!dir.exists(pr_internal)) {
    dir.create(pr_internal, recursive = TRUE)
  }

}
