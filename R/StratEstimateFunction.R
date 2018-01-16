#' Title
#'
#' @param N the geometric sampling result
#' @param cen the censored status of geometric sampling
#'
#' @return
#' @export
#'
#' @examples
StratEstimateFunction <- function(N,cen) {
  NN <- N #set d as parameters for convenience to use {boot} package later
  cend <- cen
  MLE <- rep(0,length(NN))
  MLE[cend==1] <- 1/NN[cend==1]

  result <- mean(MLE)
  return(result)
}

