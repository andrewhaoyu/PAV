
#' Title
#'
#' @param N
#' @param cen
#' @param begin
#'
#' @return
#' @export
#'
#' @examples
ParaBetaEstimateFunction <- function(N,cen,begin){
  fit <-  optim(par = begin,fn=function(x){BetaGeometricLikehood(x,N,cen)},gr=function(x){LogL.Derivatives(x,N,cen)},lower=c(0.0005,0.05),upper=c(0.9995,10000),method="L-BFGS-B",control=list(fnscale=-1))
  BetaEst <- fit$par[1]
  return(list(fit$par,BetaGeometricLikehood(fit$par,N,cen)))
}
