#' @title
#' Set up subproject's folders
#'
#' @description
#' TBD
#'
#' @param sp_name
#'
#' @return
#' @export
#'
#' @examples
setup_subproj_folders <- function(sp_name) {

  pkg_path <- get_values("pkg_path")

  sp_dir <- file.path(pkg_path, "results", sp_name)
  if (!dir.exists(sp_dir)) {dir.create(sp_dir)}

  sp_cache <- file.path(sp_dir, "cache_files")
  if (!dir.exists(sp_cache)) {dir.create(sp_cache)}

  sp_log <- file.path(sp_dir, "log_files")
  if (!dir.exists(sp_log)) {dir.create(sp_log)}

  sp_output <- file.path(sp_dir, "output_files")
  if (!dir.exists(sp_output)) {dir.create(sp_output)}

  # Add the subproject auxilary folders path into global environment.
  set_values("sp_dir", sp_dir)
  set_values("sp_cache", sp_cache)
  set_values("sp_log", sp_log)
  set_values("sp_output", sp_output)
}
