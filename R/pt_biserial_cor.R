
#' Title
#'
#' @param w TBD
#' @param binary_cov TBD
#'
#' @return
#' value
#' @export
#'
pt_biserial_cor <- function(w, binary_cov){
  # Generate point-biserial correlation between continuous exposure (w) and binary covariate (binary_cov)
  # Note - equivalent to pearson correlation
  mean_gp1 <- mean(w[binary_cov == 1])
  mean_gp0 <- mean(w[binary_cov == 0])
  s_nminus1 <- sd(w)
  n1 <- sum(binary_cov == 1)
  n0 <- sum(binary_cov == 0)
  n <- length(binary_cov) # n = n0 + n1

  cor_pb <- (mean_gp1 - mean_gp0) / s_nminus1 * sqrt(n1 / n * n0 / (n-1))
  return(cor_pb)
}
