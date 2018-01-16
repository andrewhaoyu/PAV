#' Title
#'
#' @param N the geometric sampling result
#' @param cen the censored status of geometric samplings
#'
#' @return
#' @export
#'
#' @examples
PooledEstimateFunction <- function(N,cen) {
  NN <- N
  K <- length(NN)
  cend <- cen
  result <- sum(cend)/(sum(NN)-K+sum(cend))
  return(result)
}

