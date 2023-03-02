#' @title
#' Set up required internal and external paths
#'
#' @description
#' Sets up required internal and external paths. This include any user defined
#' external data path, package path and many more.
#'
#' @param pkg A package name or path.
#'
#' @return
#' value
#' @export
#'
setup_path <- function(pkg = ".") {

  pkg <- as.package(pkg)
  pkg_path <- pkg$path

  pkg_env <- getNamespace(pkg$package)

  # collect number of valid external data path.
  val_names <- names(pkg_env)

  ext_names_pr <- grep("^DATA_PRIVATE_EXT", val_names, value = TRUE)
  ext_names_pb <- grep("^DATA_PUBLIC_EXT", val_names, value = TRUE)

  # private
  ext_names_pr_valid <- list()
  if (length(ext_names_pr) > 0){
    for (d_name in ext_names_pr){
      # get the last char after _
      tmp_str <- strsplit(d_name, "_")[[1]]
      ext_num <- tmp_str[length(tmp_str)]

      ext_num_int <- as.integer(ext_num)

      # check if it is a number, if not ignore.
      if(is.na(ext_num_int)){
        message(paste(d_name,
                      " is not a valid key for external data connection."))
        next
      }

      # check if the directory is built. If not build it.
      sub_external_name <- paste0("external_", ext_num_int)
      sub_folder_path <- file.path(pkg_path,
                                   "study_data",
                                   "private",
                                   sub_external_name)
      if (!dir.exists(sub_folder_path)) {dir.create(sub_folder_path)}
      # unlink
      # unlink(file.path(pkg_path, "study_data", "private", sub_external_name),
      #        recursive = TRUE)

      # link
      if (!is.null(pkg_env[[d_name]])){
        if (file.exists(pkg_env[[d_name]])){
          file.symlink(pkg_env[[d_name]],
                       file.path(pkg_path, "study_data", "private", sub_external_name))

          # set values.
          set_values(paste0("dir_data_private_ext_", ext_num_int),
                     file.path("study_data", "private", sub_external_name))
          ext_names_pr_valid <- c(ext_names_pr_valid,
                                  paste0("dir_data_private_ext_",
                                         ext_num_int))
        } else {
          message(paste("Provided path for '", d_name, "' is not valid."))
          next
        }
      }


    }
  }

  # public
  ext_names_pb_valid <- list()
  if (length(ext_names_pb) > 0){
    for (d_name in ext_names_pb){
      # get the last char after _
      tmp_str <- strsplit(d_name, "_")[[1]]
      ext_num <- tmp_str[length(tmp_str)]

      ext_num_int <- as.integer(ext_num)

      # check if it is a number, if not ignore.
      if(is.na(ext_num_int)){
        message(paste(d_name,
                      " is not a valid key for external data connection."))
        next
      }

      # check if the directory is built. If not build it.
      sub_external_name <- paste0("external_", ext_num_int)
      sub_folder_path <- file.path(pkg_path,
                                   "study_data",
                                   "public",
                                   sub_external_name)
      if (!dir.exists(sub_folder_path)) {dir.create(sub_folder_path)}
      # unlink
      # unlink(file.path(pkg_path, "study_data", "public", sub_external_name),
      #        recursive = TRUE)

      # link
      if (!is.null(pkg_env[[d_name]])){
        if (file.exists(pkg_env[[d_name]])){
          file.symlink(pkg_env[[d_name]],
                       file.path(pkg_path, "study_data", "public", sub_external_name))

          # set values.
          set_values(paste0("dir_data_public_ext_", ext_num_int),
                     file.path("study_data", "public", sub_external_name))
          ext_names_pb_valid <- c(ext_names_pb_valid,
                                  paste0("dir_data_public_ext_",
                                         ext_num_int))
        } else {
          message(paste("Provided path for '", d_name, "' is not valid."))
          next
        }
      }
    }
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

  set_values("ext_names_pr_valid", ext_names_pr_valid)
  set_values("ext_names_pb_valid", ext_names_pb_valid)
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
