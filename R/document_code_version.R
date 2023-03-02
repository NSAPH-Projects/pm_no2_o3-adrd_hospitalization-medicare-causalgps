#' Title
#'
#' @return
#' value
#' @export
#'
document_code_version <- function(){

 repo <- git2r::repository()

 # get the current branch name
 branch_name <- git2r::git_branch(repo)$name

 # get the current commit hash value
 commit_hash <- git2r::git_head(repo)$hash

 return(list(branch_name, commit_hash))
}
