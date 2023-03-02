#' @title
#' Initialize a subproject
#'
#' @param sp_name
#'
#' @return
#' @export
#'
#' @examples
initialize_sub_project <- function(sp_name,
                                   pkg = "."){

  # Setting up path ------------------
  setup_path(pkg = pkg)

  # Providing process name -----------
  # The process name has a three character prefix then date then any optional
  # name. We use the following convention:
  #   - inp: for in-progress exploration
  #   - man: for mature implementations that ends up providing results to include
  #          in manuscript.

  # One can review the current keys using list_key_values() function.

  setup_subproj_folders(sp_name = sp_name)

  # This concludes setting up the project. The followings are main destination to
  # use throughout the project. These variables are the same in all sub-projects
  # and should stay the same, however, it is internally set to the correct path.

  sp_dir <- get_values("sp_dir")
  sp_output <- get_values("sp_output")
  sp_cache <- get_values("sp_cache")
  sp_log <- get_values("sp_log")

  private_ext_ddir <- get_values("dir_data_private_ext")
  public_ext_ddir <- get_values("dir_data_public_ext")

  private_int_ddir <- get_values("dir_data_private_int")
  public_int_ddir <- get_values("dir_data_public_int")

  path_obj <- list()
  class(path_obj) <- "path_holder"
  path_obj$sp_dir <- sp_dir
  path_obj$sp_output <- sp_output
  path_obj$sp_cache <- sp_cache
  path_obj$sp_log <- sp_log
  path_obj$private_int_ddir <- private_int_ddir
  path_obj$public_int_ddir <- public_int_ddir

  ext_names_pr_valid <- get_values("ext_names_pr_valid")
  for (pr_path in ext_names_pr_valid){
    path_obj[pr_path] <- get_values(pr_path)
  }

  ext_names_pb_valid <- get_values("ext_names_pb_valid")
  for (pb_path in ext_names_pb_valid){
    path_obj[pb_path] <- get_values(pb_path)
  }

  return(path_obj)
}
