#' @title
#' Finalize the processing of the sub project
#'
#' @param path_obj TBD
#'
#' @return
#' value
#' @export
#'
finalize_subproject <- function(path_obj) {

  # write extra details of the session
  session_info_fn <- file.path(path_obj$sp_dir,
                               paste0("sessioninfo_",
                                      format(Sys.time(),
                                             "%Y%m%d%H%M%S"), ".txt"))
  writeLines(capture.output(sessionInfo()),
             session_info_fn)
  write("\n-------- Git:", session_info_fn, append = TRUE)
  write("\nCode version from Git:", session_info_fn, append = TRUE)
  write("\nCurrent HEAD info:", session_info_fn, append = TRUE)
  write(capture.output(git2r::repository_head()), session_info_fn, append = TRUE)
  write("\nLast commit info:", session_info_fn, append = TRUE)
  for (i in seq(1,2)){
    write(paste(git2r::commits()[[1]][i]), session_info_fn, append = TRUE)
  }

  write("\n-------- Path:", session_info_fn, append = TRUE)
  write("\nContent of R/external_path.R:", session_info_fn, append = TRUE)
  write(readLines("R/external_path.R"), session_info_fn, append = TRUE)

  write("\nContent of path_obj:", session_info_fn, append = TRUE)
  write(capture.output(print(path_obj)), session_info_fn, append = TRUE)

}
