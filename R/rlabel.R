#' Generation of a missing-data indicator
#'
#' Generate the missing label indicator
#'
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation.
#' @param pi A g-dimensional  initial vector of the mixing proportions.
#' @param mu A initial  \eqn{p \times g} matrix of the location parameters.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @param xi A 2-dimensional coefficient vector for a logistic function of the Shannon entropy.
#' @return
#' \item{m}{A n-dimensional vector of missing label indicator. The element of  outputs \code{m} represents its label indicator is missing if m equals 1, otherwise its label indicator is available if m equals to 0.}
#' @export
#' @import stats
#' @examples
#' n<-150
#' pi<-c(0.25,0.25,0.25,0.25)
#' sigma<-array(0,dim=c(3,3,4))
#' sigma[,,1]<-diag(1,3)
#' sigma[,,2]<-diag(2,3)
#' sigma[,,3]<-diag(3,3)
#' sigma[,,4]<-diag(4,3)
#' mu<-matrix(c(0.2,0.3,0.4,0.2,0.7,0.6,0.1,0.7,1.6,0.2,1.7,0.6),3,4)
#' dat<-rmix(n=n,pi=pi,mu=mu,sigma=sigma)
#' xi<-c(-0.5,1)
#' m<-rlabel(dat=dat$Y,pi=pi,mu=mu,sigma=sigma,xi=xi)

rlabel <- function(dat, pi, mu, sigma,xi){
  X=dat
  n <- nrow(as.matrix(X))
  g<-length(pi)
  p<-dim(mu)[1]
  ncov=ifelse(is.na(dim(sigma)[3]),1,2)
  if(ncov==1){
    if(g==2){
      betapar <- discriminant_beta(pi, mu, sigma)
      beta0 <- betapar$beta0
      beta <- betapar$beta
      omega <- beta0 + X %*% beta  # should be t(beta)%*%X
      eta <- xi[1] + xi[2]*omega^2  #replace the entropy by the square the discriminant function d()
      rmu <- 1/(1+exp(-eta))
      m <- stats::rbinom(n, size=1, prob=rmu)
    }else{
      sigma1<-array(0,dim=c(p,p,g))
      for(i in 1:g){
        sigma1[,,i]=sigma
      }
      dfun=get_entropy(dat=X, n=n, p=p, g=g,mu=mu, sigma=sigma1,pi=pi)
      eta=xi[1]+xi[2]*(log(dfun))
      rmu <- 1/(1+exp(-eta))
      m <- stats::rbinom(n, size=1, prob=rmu)
    }
  }else{
    dfun=get_entropy(dat=X, n=n, p=p, g=g, mu=mu, sigma=sigma,pi=pi)
    eta=xi[1]+xi[2]*(log(dfun))
    rmu <- 1/(1+exp(-eta))
    m <- stats::rbinom(n, size=1, prob=rmu)
  }
  return(m)
}
