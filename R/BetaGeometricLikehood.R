#' Title
#'
#' @param x
#' @param N
#' @param cen
#'
#' @return
#' @export
#'
#' @examples
BetaGeometricLikehood <- function(x,N,cen){

  ThetaBar.Par<- x[1]
  M.Par <- x[2]
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- lgamma(alpha+beta)+lgamma(alpha+cen)+lgamma(beta+N-1)-lgamma(alpha)-lgamma(beta)-lgamma(alpha+beta+N+cen-1)
  return(sum(result))
}



#' Title
#'
#' @param x
#' @param N
#' @param cen
#'
#' @return
#' @export
#'
#' @examples
LogL.Derivatives <- function(x,N,cen){
  return(c(LogL.Dmu(x[1],x[2],N,cen),LogL.DM(x[1],x[2],N,cen)))
}

#' Title
#'
#' @param ThetaBar.Par
#' @param M.Par
#' @param N
#' @param cen
#'
#' @return
#' @export
#'
#' @examples
LogL.Dmu <- function(ThetaBar.Par,M.Par,N,cen){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- LogL.Dalpha(ThetaBar.Par,M.Par,N,cen)*M.Par-M.Par*LogL.Dbeta(ThetaBar.Par,M.Par,N,cen)
  return(result)
}


#' Title
#'
#' @param ThetaBar.Par
#' @param M.Par
#' @param N
#' @param cen
#'
#' @return
#' @export
#'
#' @examples
LogL.DM <- function(ThetaBar.Par,M.Par,N,cen){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- LogL.Dalpha(ThetaBar.Par,M.Par,N,cen)*ThetaBar.Par-(1-ThetaBar.Par)*LogL.Dbeta(ThetaBar.Par,M.Par,N,cen)
}



#' Title
#'
#' @param ThetaBar.Par
#' @param M.Par
#' @param N
#' @param cen
#'
#' @return
#' @export
#'
#' @examples
LogL.Dalpha <- function(ThetaBar.Par,M.Par,N,cen){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- digamma(alpha+beta)/(exp(lgamma(alpha+beta)))+digamma(alpha+cen)/exp(lgamma(alpha+cen))-digamma(alpha)/exp(lgamma(alpha))
  -digamma(alpha+beta+N+cen-1)/exp(lgamma(alpha+beta+N+cen-1))
  return(sum(result))
}


#' Title
#'
#' @param ThetaBar.Par
#' @param M.Par
#' @param N
#' @param cen
#'
#' @return
#' @export
#'
#' @examples
LogL.Dbeta <- function(ThetaBar.Par,M.Par,N,cen){
  alpha <- ThetaBar.Par*M.Par
  beta <- M.Par*(1-ThetaBar.Par)
  result <- digamma(alpha+beta)/exp(lgamma(alpha+beta))+digamma(beta+N-1)/exp(lgamma(beta+N-1))-
    digamma(beta)/exp(lgamma(beta))-digamma(alpha+beta+N+cen-1)/exp(lgamma(alpha+beta+N+cen-1))
  sum(result)
}
