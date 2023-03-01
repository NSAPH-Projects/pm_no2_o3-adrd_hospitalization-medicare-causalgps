#' Title
#'
#' @param df
#'
#' @return
#' @export
#'
#' @examples
explore_zip_covs <- function(df){
  # Explore distribution of ZIP-level covariates
  cat("\nMin of mean_bmi:", min(df$mean_bmi))
  cat("\nMax of mean_bmi:", max(df$mean_bmi))
  cat("\nMin of smoke_rate:", min(df$smoke_rate))
  cat("\nMax of smoke_rate:", max(df$smoke_rate))
  cat("\nMin of prop_blk:", min(df$prop_blk))
  cat("\nMax of prop_blk:", max(df$prop_blk))
  cat("\nMin of hispanic:", min(df$hispanic))
  cat("\nMax of hispanic:", max(df$hispanic))
  cat("\nMin of education:", min(df$education))
  cat("\nMax of education:", max(df$education))
  cat("\nMin of popdensity:", min(df$popdensity))
  cat("\nMax of popdensity:", max(df$popdensity))
  cat("\nMin of poverty:", min(df$poverty))
  cat("\nMax of poverty:", max(df$poverty))
  cat("\nMin of PIR:", min(df$PIR))
  cat("\nMax of PIR:", max(df$PIR))
  cat("\nMin of prop_owner_occ:", min(df$prop_owner_occ))
  cat("\nMax of prop_owner_occ:", max(df$prop_owner_occ))
  cat("\nMin of summer_tmmx:", min(df$summer_tmmx))
  cat("\nMax of summer_tmmx:", max(df$summer_tmmx))
  cat("\nMin of summer_rmax:", min(df$summer_rmax))
  cat("\nMax of summer_rmax:", max(df$summer_rmax))
  cat("\nMin of no2:", min(df$no2))
  cat("\nMax of no2:", max(df$no2))
  cat("\nMin of ozone_summer:", min(df$ozone_summer))
  cat("\nMax of ozone_summer:", max(df$ozone_summer))
  prop.table(table(df$region))

  # cat("\nMean of mean_bmi:", mean(df$mean_bmi))
  # cat("\nSD of mean_bmi:", sd(df$mean_bmi))
  # cat("\nMean of smoke_rate:", mean(df$smoke_rate))
  # cat("\nSD of smoke_rate:", sd(df$smoke_rate))
  # cat("\nMean of prop_blk:", mean(df$prop_blk))
  # cat("\nSD of prop_blk:", sd(df$prop_blk))
  # cat("\nMean of hispanic:", mean(df$hispanic))
  # cat("\nSD of hispanic:", sd(df$hispanic))
  # cat("\nMean of education:", mean(df$education))
  # cat("\nSD of education:", sd(df$education))
  # cat("\nMean of popdensity:", mean(df$popdensity))
  # cat("\nSD of popdensity:", sd(df$popdensity))
  # cat("\nMean of poverty:", mean(df$poverty))
  # cat("\nSD of poverty:", sd(df$poverty))
  # cat("\nMean of medhouseholdincome:", mean(df$medhouseholdincome))
  # cat("\nSD of medhouseholdincome:", sd(df$medhouseholdincome))
  # cat("\nMean of PIR:", mean(df$PIR))
  # cat("\nSD of PIR:", sd(df$PIR))
  # cat("\nMean of prop_owner_occ:", mean(df$prop_owner_occ))
  # cat("\nSD of prop_owner_occ:", sd(df$prop_owner_occ))
}
