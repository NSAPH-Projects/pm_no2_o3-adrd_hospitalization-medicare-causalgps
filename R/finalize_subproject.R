#' @title
#' Finalize the processing of the sub project
#'
#' @param path_obj
#'
#' @return
#' @export
#'
#' @examples
finalize_subproject <- function(path_obj) {

  # write extra details of the session
  session_info_fn <- file.path(path_obj$sp_output,
                               paste0("sessioninfo_",
                                      format(Sys.time(),
                                             "%Y%m%d%H%M%S"), ".txt"))
  writeLines(capture.output(sessionInfo()),
             session_info_fn)
  write("\nCode version from Git:", session_info_fn, append = TRUE)
  for (i in seq(1,2)){
    write(paste(git2r::commits()[[1]][i]), session_info_fn, append = TRUE)
  }

}
