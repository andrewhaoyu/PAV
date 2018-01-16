#' Title
#'
#' @param N the geometric sampling result
#'
#' @return
#' @export
#'
#' @examples
ParaBetaEstimateFunction <- function(N){
  fit <-  optim(par = begin,fn=BetaGeometricLikehood,gr=LogL.Derivatives,lower=c(0.0005,0.05),upper=c(0.5,10000),method="L-BFGS-B",control=list(fnscale=-1))

  return(fit$par[1])
}
