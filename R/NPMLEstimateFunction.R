
#' Title
#'
#' @param N: the geometric sampling result
#' @param cen: the censored status of geometric sampling
#' @param d: an index to select part of the data, default is 1:N. make it easy for the boot package.
#'
#' @return the discrete probability mass of the underlying distribution. UU is the probability mass, pp is the corresponding probability on that probability mass
#' @export
#'
#' @examples
NPMLEstimateFunction <- function(N,cen,d=c(1:length(N))) {
  #beginning setting of NPML:
  # NN is the geometric sampling result
  # UU is the probability mass point of NPML
  # UU begins with 1:K/(K+1)
  # ww is the posterior probability distribution of ith sample coming from UU_jth
  # YY is the a vector with all elements equals 1,
  # YY setting is easy for negative bionomial case
  # pp is the probability distrubution of UU
  #recursion steps of NPML:
  # 500 Niter steps for the EM algorithm
  # a standard to break the loop,
  # if sum(UU_nextstep-UU_now)^2<1e-5&sum(ww_nextstep-ww_now)^2<1e-5&steps>25
  # set a epi standard for UU, if UU>(1-epi),UU=1-epi; if UU<epi, UU=epi;
  #results of the NPML:
  # we return a lis of UU and pp for the NPML
  NN <- N[d]
  KK <- length(NN)
  cend = cen[d]
  Y <- rep(1,KK)
  UU <- 1:KK/(KK+1)

  ww <- matrix(1, nrow = KK, ncol=KK)
  pp <- rep(1/KK, length = KK)
  Niter <- 500
  temp_UU=1
  epi=0.001
  temp_pp=1
  tol <- 1e-5
  # The recursion

  for(i in 1:Niter) {
    #Niter

    for(k in 1:KK) {
      ww[k,] <- mapply(function(x,y) recursion(x,y,NN[k],Y[k],cend[k]),pp,UU)
      ww[k,] <- ww[k,]/sum(ww[k,])
    }  #k to K
    pp <- apply(ww,2,sum)/sum(ww)
    for(j in 1:KK) {
      UU[j] <- sum(ww[,j]*Y*cend)/sum(ww[,j]*(NN-Y+cend))
    }
    UU[which(UU>(1-epi))] <- 1-epi
    UU[which(UU<epi)] <- epi
    if(any(is.nan(UU))==T) {
      break
    }
    if((sum((temp_UU-UU)^2)<tol)&sum((temp_pp-pp)^2)<tol) {
      break
    }
    temp_pp <- pp
    temp_UU <- UU
  }


  return(list(UU,pp,sum(UU*pp)))

}


