
setup_path <- function(pkg = "."){

  pkg <- as.package(pkg)
  pkg_path <- pkg$path

  # Load from .Renviron file
  readRenviron(".Renviron")

  # Collect values for keys
  dir_d_private_ext <- Sys.getenv("DATA_PRIVATE_EXT")
  dir_d_public_ext <- Sys.getenv("DATA_PUBLIC_EXT")

  dir_d_private_ext_name <- basename(dir_d_private_ext)
  dir_d_public_ext_name <- basename(dir_d_public_ext)

  # Unlink possible available links
  unlink(file.path(pkg_path, "data", "public", "external"), recursive = TRUE)
  unlink(file.path(pkg_path, "data", "private", "external"), recursive = TRUE)

  file.symlink(dir_d_public_ext, file.path(pkg_path,
                                             "data", "public", "external"))
  file.symlink(dir_d_private_ext, file.path(pkg_path,
                                               "data", "private", "external"))


  # Set up values
  set_values("dir_data_private_ext", file.path("data", "private", "external",
                                               dir_d_private_ext_name))
  set_values("dir_data_public_ext", file.path("data", "public", "external",
                                               dir_d_public_ext_name))

  set_values("dir_data_private_internal", file.path(pkg_path,
                                                    "data",
                                                    "private",
                                                    "internal"))

  set_values("dir_data_public_internal", file.path(pkg_path,
                                                   "data",
                                                   "public",
                                                   "internal"))

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
