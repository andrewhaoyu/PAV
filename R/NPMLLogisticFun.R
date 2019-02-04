#' logisitic function
#'
#' @param x
#'
#' @return
logit_inver <- function(x){
  return(exp(x)/(1+exp(x)))
}

#' E step for NPML logisitic
#'
#' @param uu
#' @param beta
#' @param x
#' @param y
#' @param w
#'
#' @return
Estep <- function(uu,beta,x,y,w){
  n <- length(y)
  ###let i represent person
  ###let j represent the population unit
  ###pp calculate the pr(N_i|Z_j) = (1-p_j)^(N_i-1)*(1-p_j)

  pp <- matrix(0,n,n)
  xb <- x%*%beta
  for(i in 1:K){
    for(j in 1:K){
      pp[i,j] <- (1-logit_inver(uu[j]+xb[i]))^(y[i]-1)*logit_inver(uu[j]+xb[i])
    }
  }
  ######### ww perform the bayesian rule
  ######### Pr(Z_j|N_i) = Pr(N_i|Z_j)*Pr(Z_j)/(sum_j Pr(N_i|Z_j)*Pr(Z_j))

  ww <- matrix(0,n,n)
  for(i in 1:K){

    for(j in 1:K){
      #ww right now is the joint distribution of i,j
      ww[i,j] = w[j]*pp[i,j]
    }
  }
  ##posterior = joint distribution/marginal distribution
  for(i in 1:K){
    ww[i,] = ww[i,]/sum(ww[i,])
  }

  return(ww)
}

#' gradient decent for M step intercept
#'
#' @param u
#' @param N
#' @param ww
#' @param b
#' @param x
#' @param K
#'
#' @return
gr_u_fun <- function(u,N,ww,b,x,K){
  sum <- 0
  xb <- x%*%b
  gr_u <- rep(0,K)
  #record a temp matrix (1-n_k*p)
  temp_mat <- matrix(0,K,K)
  for(i in 1:K){
    temp_mat[i,] = 1-N[i]*logit_inver(u+xb[i])

  }
  #use the element wise times in R
  gr_u = colSums(temp_mat*ww)
  # for(j in 1:K){
  #   gr_u[j] <- crossprod(temp_mat[,j],ww[,j])
  #   }
  return(gr_u)
}

#' gradient decent for M step slope
#'
#' @param uu
#' @param N
#' @param ww
#' @param b
#' @param x
#' @param K
#'
#' @return
#' @export
#'
#' @examples
gr_b_fun <- function(uu,N,ww,b,x,K){
  u <- uu
  p <- length(b)
  t <- 0
  xb <- x%*%b
  temp_mat <- matrix(0,K,K)
  for(i in 1:K){
    temp_mat[i,] = 1-N[i]*logit_inver(u+xb[i])

  }

  for(k in 1:p){
    x_temp <- x[,k]
    #use the element times in R
    p[k] <- sum(x_temp*temp_mat*ww)
  }

  return(p)
}


#' M Step for NPMLLogfun
#'
#' @param uu
#' @param beta
#' @param x
#' @param y
#' @param ww
#'
#' @return
Mstep <- function(uu,beta,x,y,ww){
  uu_old <- uu
  beta_old <- beta
  step <- 100
  n <- length(y)
  alpha <- 1/n
  tol <- 0.001
  for(l in 1:step){

    uu_beta_old <- c(uu_old,beta_old)
  #print(uu_beta_old)
    uu_new = uu_old+alpha*gr_u_fun(uu_old,y,ww,beta_old,x,n)
    beta_new <- beta_old+alpha*gr_b_fun(uu_new,y,ww,beta_old,x,n)
    uu_beta_new <- c(uu_new,beta_new)
    error <- max(abs(uu_beta_new-uu_beta_old))
    if(error<tol){
      break
    }
    uu_old <- uu_new
    beta_old <- beta_new
  }
  return(list(uu_new,beta_new))
}


#' NPMLLogFun
#'
#' @param y
#' @param x
#' @param uu0
#' @param beta0
#'
#' @return
#' @export
#'
NPMLLogFun <- function(y,x,uu0,beta0){
  x <- as.matrix(x)
  step = 100
  uu_old = uu0
  beta_old = beta0
  tol = 0.001
  n = length(y)
  w = rep(1/n,n)
  for(l in 1:step){
    uu_beta_old <- c(uu_old,beta_old)
    print(uu_beta_old)
    ww = Estep(uu_old,beta_old,x,y,w)
    Mstep_result = Mstep(uu_old,beta_old,x,y,ww)
    uu_new = Mstep_result[[1]]
    beta_new = Mstep_result[[2]]
    uu_beta_new <- c(uu_new,beta_new)
    error <- max(abs(uu_beta_new-uu_beta_old))
    if(error<tol){
      break
    }
    uu_old <- uu_new
    beta_old <- beta_new
    w = colSums(ww)/sum(ww)
  }
  return(list(uu_new,beta_new))
}
