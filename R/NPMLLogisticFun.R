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
#' @param obs
#'
#' @return
Estep <- function(uu,beta,x,y,w,obs){
  n <- length(y)
  ###let i represent person
  ###let j represent the population unit
  ###pp calculate the pr(N_i|Z_j) = (1-p_j)^(N_i-1)*(1-p_j)

  pp <- matrix(0,n,n)
  xb <- x%*%beta
  for(i in 1:n){
    p_temp = logit_inver(uu+xb[i])
pp[i,] <- (1-p_temp)^(y[i]-obs[i])*p_temp^obs[i]
}
  ######### ww perform the bayesian rule
  ######### Pr(Z_j|N_i) = Pr(N_i|Z_j)*Pr(Z_j)/(sum_j Pr(N_i|Z_j)*Pr(Z_j))

  ww <- matrix(0,n,n)
  for(i in 1:n){
#ww right now is the joint distribution of i,j
    ww[i,] = w*pp[i,]
  }
  ##posterior = joint distribution/marginal distribution
  for(i in 1:n){
    ww[i,] = ww[i,]/sum(ww[i,])
  }

  return(ww)
}




#' Title
#'
#' @param y
#' @param x
#' @param uu_old
#' @param beta_old
#' @param w
#' @param obs
#'
#' @return
ObsLikfun <- function(y,x,uu_old,beta_old,w,obs){
  n <- length(y)
  result <- rep(0,n)
  xb = x%*%beta_old
  for(i in 1:n){
    result[i] <- log(sum(w*logit_inver(uu_old+xb[i])^obs[i]*(1-logit_inver(uu_old+xb[i]))^(y[i]-obs[i])))
  }
  return(sum(result))

}

#' Title
#'
#' @param y
#' @param x
#' @param uu_old
#' @param beta_old
#' @param ww
#' @param obs
#'
#' @return
ComLikfun <- function(y,x,uu_old,beta_old,ww,obs){
  xb = x%*%beta_old
  n <- length(y)
  temp_mat = matrix(0,n,n)
  for(i in 1:n){
    temp_mat[i,] = ((uu_old+xb[i])*obs[i]-y[i]*log(1+exp(uu_old+xb[i])))
  }
  return(sum(temp_mat*ww))
}





#' gradient decent for M step intercept
#'
#' @param u
#' @param N
#' @param ww
#' @param b
#' @param x
#' @param K
#' @param obs
#'
#' @return
gr_u_fun <- function(u,N,ww,b,x,K,obs){
  sum <- 0
  xb <- x%*%b
  gr_u <- rep(0,K)
  #record a temp matrix (1-n_k*p)
  temp_mat <- matrix(0,K,K)
  for(i in 1:K){
    temp_mat[i,] = obs[i]-N[i]*logit_inver(u+xb[i])
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
#' @param obs
#'
#' @return
#' @export
#'
#' @examples
gr_b_fun <- function(uu,N,ww,b,x,K,obs){
  u <- uu
  p <- length(b)
  t <- 0
  xb <- x%*%b
  temp_mat <- matrix(0,K,K)
  for(i in 1:K){
    temp_mat[i,] = obs[i]-N[i]*logit_inver(u+xb[i])

  }

  for(k in 1:p){
    x_temp <- x[,k]
    #use the element times in R
    p[k] <- sum(x_temp%*%(temp_mat*ww))
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
#' @param alpha_x
#' @param obs
#'
#' @return
Mstep <- function(uu_old,beta_old,x,y,ww,alpha_x,obs){
  step <- 200
  n <- length(y)
  alpha <- 1/(n*10)
  tol <- alpha
  ComLik0 = ComLikfun(y,x,uu_old,beta_old,ww,obs)
  #Likelihood_old <- Likelihoodfun(y,x,uu_old,beta_old,ww)
  for(l in 1:step){

    uu_beta_old <- c(uu_old,beta_old)
    #print(uu_beta_old)
    uu_new = uu_old+alpha*gr_u_fun(uu_old,y,ww,beta_old,x,n,obs)
    beta_new <- beta_old+(alpha_x/100)*gr_b_fun(uu_new,y,ww,beta_old,x,n,obs)
    ######test the Likelihood for a few steps
    if(l%%10==0){
      ComLik_new = ComLikfun(y,x,uu_new,beta_new,ww,obs)
      if(ComLik_new>ComLik0){
        break
      }else if(ComLik_new<=ComLik0){
        uu_new = uu_old
        beta_new = beta_old
        alpha = alpha/2
        alpha_x = alpha_x/2
      }
    }
    uu_beta_new <- c(uu_new,beta_new)
    error <- max(abs(uu_beta_new-uu_beta_old))
    if(error<tol){
      break
    }
    uu_old <- uu_new
    beta_old <- beta_new
  }
  return(list(uu_new,beta_new,l))
}



#' NPMLLogFun
#'
#' @param y
#' @param x
#' @param obs
#' @param uu0
#' @param beta0
#'
#' @return
#' @export
#'
NPMLLogFun <- function(y,x,obs,uu0,beta0){
  if(is.null(uu0)){
    #use the stratified propobablity with some smoother
    y_sm = y+1
    uu0=seq(min(log((1/y_sm)/(1-1/y_sm))),
            max(log((1/y_sm)/(1-1/y_sm))),
            max(log((1/y_sm)/(1-1/y_sm)))-min(log((1/y_sm)/(1-1/y_sm)))/(n-1))
  }
  x <- as.matrix(x)
  step = 100
  uu_old = uu0
  beta_old = beta0
  tol = 0.001
  n = length(y)
  w = rep(1/n,n)
  #set the step length of gradient decent to aviod unconvergence
  var.x = apply(x,2,var)
  alpha_x = rep(1/n,length(beta_old))
  for(i in 1:length(beta_old)){
  if(var.x[i]>=1){
    alpha_x[i] = alpha_x[i]/var.x[i]
  }
  }
  LikeliResult <- rep(0,step)
  StepResult <- rep(0,step)
  #library(PAV)
  for(l in 1:step){
    uu_beta_old <- c(uu_old,beta_old)
    #print(uu_beta_old)
    #print(uu_beta_old)
    ww = Estep(uu_old,beta_old,x,y,w,obs)
    LikeliResult[l] <- ObsLikfun(y,x,uu_old,beta_old,w,obs)
    #rowSums(ww)
    Mstep_result = Mstep(uu_old,beta_old,x,y,ww,alpha_x,obs)
    #Mstep_result = Mstep2(uu_old,beta_old,x,y,ww)
    uu_new = Mstep_result[[1]]
    beta_new = Mstep_result[[2]]
    #StepResult[l] = Mstep_result[[3]]
    uu_beta_new <- c(uu_new,beta_new)
    error <- max(abs(uu_beta_new-uu_beta_old))
    if(error<tol){
      break
    }
    uu_old <- uu_new
    beta_old <- beta_new
    w = colSums(ww)/sum(ww)
  }
  max_likelihood <-   LikeliResult[l]
  max_step = l

  return(list(uu_new,beta_new,max_likelihood,max_step))
}
