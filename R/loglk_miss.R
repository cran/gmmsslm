#' Log likelihood function formed on the basis of the missing-label indicator
#'
#' Log likelihood for partially classified data based on the missing mechanism with the Shanon entropy
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param zm An n-dimensional vector containing the class labels including the missing-label denoted as NA.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix if \code{ncov=1}, or a list of g covariance matrices with dimension \eqn{p\times p \times g} if \code{ncov=2}.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @param xi A 2-dimensional vector containing the initial values of the coefficients in the logistic function of the Shannon entropy.
#' @return
#' \item{lk}{loglikelihood value}
#' @export
#'
#' @details The log-likelihood function  formed on the basis of the missing-label indicator can be expressed by
#' \deqn{
#' \log L_{PC}^{({miss})}(\theta,\boldsymbol{\xi})=\sum_{j=1}^n\big[ (1-m_j)\log\left\lbrace 1-q(y_j;\theta,\boldsymbol{\xi})\right\rbrace +m_j\log q(y_j;\theta,\boldsymbol{\xi})\big],
#' }
#'where \eqn{q(y_j;\theta,\boldsymbol{\xi})} is a logistic function of the Shannon entropy \eqn{e_j(y_j;\theta)},
#'and  \eqn{m_j} is a missing label indicator.
#'
#'
loglk_miss<- function(dat,zm,pi,mu,sigma,xi){
  Y<-dat
  n <- nrow(Y)
  g<-length(pi)
  p<-dim(mu)[1]
  ncov=ifelse(is.na(dim(sigma)[3]),1,2)
  if(ncov==1){
    if(g==2){
      betapar <- discriminant_beta(pi, mu, sigma)
      beta0 <- betapar$beta0
      beta <- betapar$beta
      omega <- beta0 + as.vector(Y %*% beta)
      eta <- xi[1] + xi[2]*omega^2
    }else{
      sigma1<-array(0,dim=c(p,p,g))
      for(i in 1:g){
        sigma1[,,i]=sigma
      }
      dfun=get_entropy(dat=Y, n=n, p=p, g=g,  mu=mu, sigma=sigma1,pi=pi)
      eta=xi[1]+xi[2]*(log(dfun))
    }
  }else{
    dfun=get_entropy(dat=Y, n=n, p=p, g=g, mu=mu, sigma=sigma,pi=pi)
    eta=xi[1]+xi[2]*(log(dfun))   ###############change the square of the discriminant funtion to to entropy
    #    }
  }
  m <- as.numeric((is.na(zm)))
  n1 <- sum(m==1)
  n2 <- sum(m==0)
  if (n1 > 0){
    lk1 <- sum(-sapply(eta[m==1], function(x) logsumexp(c(0, -x))))
  } else {
    lk1 <- 0
  }
  if (n2 > 0){
    lk0 <- sum(-sapply(eta[m==0], function(x) logsumexp(c(0, x))))
  } else {
    lk0 <- 0
  }
  lk <- lk1+lk0
  return(lk)
}
