#' Title
#'
#' @param x the geometric sampling result
#'
#' @return
#' @export
#'
#' @examples
BetaGeometricLikehood <- function(x){

  ThetaBar.Par<- x[1]
  M.Par <- x[2]
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- lgamma(alpha+beta)+lgamma(alpha+cen)+lgamma(beta+N-1)-lgamma(alpha)-lgamma(beta)-lgamma(alpha+beta+N+cen-1)
  return(sum(result))
}



#' Title
#'
#' @param x the geometric sampling result
#'
#' @return
#' @export
#'
#' @examples
LogL.Derivatives <- function(x){
  return(c(LogL.Dmu(x[1],x[2]),LogL.DM(x[1],x[2])))
}

#' Title
#'
#' @param ThetaBar.Par Beta distribution mean paramter
#' @param M.Par Beta distribution precision paramter
#'
#' @return
#' @export
#'
#' @examples
LogL.Dmu <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- LogL.Dalpha(ThetaBar.Par,M.Par)*M.Par-M*LogL.Dbeta(ThetaBar.Par,M.Par)
  return(result)
}

#' Title
#'
#' @param ThetaBar.Par Beta distribution mean paramter
#' @param M.Par Beta distribution precision paramter
#'
#' @return
#' @export
#'
#' @examples
LogL.DM <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- LogL.Dalpha(ThetaBar.Par,M.Par)*ThetaBar.Par-(1-ThetaBar.Par)*LogL.Dbeta(ThetaBar.Par,M.Par)
}

#' Title
#'
#' @param ThetaBar.Par Beta distribution mean paramter
#' @param M.Par Beta distribution precision paramter
#'
#' @return
#' @export
#'
#' @examples
LogL.Dalpha <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- digamma(alpha+beta)/(exp(lgamma(alpha+beta)))+digamma(alpha+cen)/exp(lgamma(alpha+cen))-digamma(alpha)/exp(lgamma(alpha))
  -digamma(alpha+beta+N+cen-1)/exp(lgamma(alpha+beta+N+cen-1))
  return(sum(result))
}

#' Title
#'
#' @param ThetaBar.Par Beta distribution mean paramter
#' @param M.Par Beta distribution precision paramter
#'
#' @return
#' @export
#'
#' @examples
LogL.Dbeta <- function(ThetaBar.Par,M.Par){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- digamma(alpha+beta)/exp(lgamma(alpha+beta))+digamma(beta+N-1)/exp(lgamma(beta+N-1))-
    digamma(beta)/exp(lgamma(beta))-digamma(alpha+beta+N+cen-1)/exp(lgamma(alpha+beta+N+cen-1))
  sum(result)
}
