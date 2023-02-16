#'  Shannon entropy
#'
#'  Shannon entropy
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param n Number of observations.
#' @param p Dimension of observation vecor.
#' @param g Number of multivariate normal classes.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @return
#' \item{clusprobs}{The posterior probabilities of the i-th entity that belongs to the j-th group.}
#' @export
#' @import mvtnorm
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
#' en<-get_entropy(dat=dat$Y,n=150,p=3,g=4,mu=mu,sigma=sigma,pi=pi)
#' @details
#' The concept of information entropy was introduced by \cite{shannon1948mathematical}.
#' The entropy of \eqn{y_j} is formally defined as
#' \deqn{e_j( y_j; \theta)=-\sum_{i=1}^g \tau_i( y_j; \theta) \log\tau_i(y_j;\theta).}




get_entropy <- function(dat, n, p, g, pi, mu, sigma){
  tau <- get_clusterprobs(dat=dat, n=n, p=p, g=g, mu=mu, sigma=sigma,pi=pi)
  entropy <- apply(tau, 1, function(x) sum(ifelse(x==0, 0,-x*log(x))))
  return(entropy)
}
