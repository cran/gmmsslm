#' Log likelihood for partially classified data with ingoring the missing mechanism
#'
#' Log likelihood for partially classified data with ingoring the missing mechanism
#'
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param zm An n-dimensional vector containing the class labels including the missing-label denoted as NA.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#'
#' @return
#'  \item{lk}{Log-likelihood value.}
#' @details
#'  The log-likelihood function for  partially classified data with ingoring the missing mechanism can be expressed as
#'  \deqn{
#'  \log L_{PC}^{({ig})}(\theta)=\sum_{j=1}^n  \left[
#' (1-m_j)\sum_{i=1}^g z_{ij}\left\lbrace \log\pi_i+\log  f_i(y_j;\omega_i)\right\rbrace +m_j\log \left\lbrace  \sum_{i=1}^g\pi_i  f_i(y_j;\omega_i)\right\rbrace  \right],
#'  }
#'  where \eqn{m_j} is a missing label indicator, \eqn{z_{ij}} is a zero-one indicator variable defining the known group of origin of each,
#'  and \eqn{f_i(y_j;\omega_i)} is a probability density function with parameters \eqn{\omega_i}.
#'
#'
#'
#' @export

loglk_ig <- function(dat,zm,pi, mu, sigma){
  Y<-dat
  ncov=ifelse(is.na(dim(sigma)[3]),1,2)
  g<-length(pi)
  ni<-NULL
  grp<-NULL
  for(i in 1:g){
    ni[i]<-sum(zm==i, na.rm=TRUE)
  }
  nmiss <- sum(is.na(zm))
  if(ncov==1){
    for(ii in 1:g){
      if(ni[i]>0){
        grp[ii] <- sum(mvtnorm::dmvnorm(Y[which(zm==ii),,drop=FALSE], mean=mu[,ii,drop=FALSE], sigma=as.matrix(sigma), log=TRUE))+ni[ii]*log(pi[ii])
      }else{
        grp[ii]<-0
      }
    }
    lablk<-sum(grp)
    if(nmiss>0){
      D <- matrix(0, nmiss,g)
      for(j in 1:g){
        D[,j] <- mvtnorm::dmvnorm(Y[which(is.na(zm)),,drop=FALSE], mean=mu[,j,drop=FALSE], sigma=as.matrix(sigma), log=TRUE)+log(pi[j])
      }
      unlablk <- sum(apply(D, 1, logsumexp))
    }else {
      unlablk <- 0
    }
  }else{
    for(ii in 1:g){
      if(ni[i]>0){
        grp[ii] <- sum(mvtnorm::dmvnorm(Y[which(zm==ii),,drop=FALSE], mean=mu[,ii,drop=FALSE], sigma=as.matrix(sigma[,,ii]), log=TRUE))+ni[ii]*log(pi[ii])
      }else{
        grp[ii]<-0
      }
    }
    lablk<-sum(grp)
    if(nmiss>0){
      D <- matrix(0, nmiss,g)
      for(j in 1:g){
        D[,j] <- mvtnorm::dmvnorm(Y[which(is.na(zm)),,drop=FALSE], mean=mu[,j,drop=FALSE], sigma=as.matrix(sigma[,,j]), log=TRUE)+log(pi[j])
      }
      unlablk <- sum(apply(D, 1, logsumexp))
    }else {
      unlablk <- 0
    }
  }
  lk <- lablk+unlablk

  return(lk)
}
