#' @title
#' Set up required internal and external paths
#'
#' @description
#' Sets up required internal and external paths. This include any user defined
#' external data path, package path and many more.
#'
#' @param pkg
#' @param external_private
#' @param external_public
#'
#' @return
#' @export
#'
#' @examples
setup_path <- function(pkg = ".",
                       external_private = TRUE,
                       external_public = TRUE) {

  pkg <- as.package(pkg)
  pkg_path <- pkg$path

  if (external_private){
    dir_d_private_ext <- DATA_PRIVATE_EXT
    unlink(file.path(pkg_path, "study_data", "private", "external"),
           recursive = TRUE)
    file.symlink(dir_d_private_ext,
                 file.path(pkg_path, "study_data", "private", "external"))
    set_values("dir_data_private_ext",
               file.path("study_data", "private", "external"))
  }

  if (external_public){
    dir_d_public_ext <- DATA_PUBLIC_EXT
    unlink(file.path(pkg_path, "study_data", "public", "external"),
           recursive = TRUE)
    file.symlink(dir_d_public_ext,
                 file.path(pkg_path, "study_data", "public", "external"))
    set_values("dir_data_public_ext",
               file.path("study_data", "public", "external"))
  }

  set_values("dir_data_private_int", file.path(pkg_path,
                                                    "study_data",
                                                    "private",
                                                    "internal"))

  set_values("dir_data_public_int", file.path(pkg_path,
                                                   "study_data",
                                                   "public",
                                                   "internal"))

  set_values("pkg_path", pkg_path)
}

internal_params <- new.env(parent = emptyenv())

get_values <- function(key) {
  internal_params[[key]]
}

set_values <- function(key, value) {
  internal_params[[key]] <- value
}

list_key_values <- function() {
  names(internal_params)
}
