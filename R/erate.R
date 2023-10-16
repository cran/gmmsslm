#' Error rate of the Bayes rule for a g-class Gaussian mixture model
#'
#'Error rate of the Bayes rule for a g-class Gaussian mixture model
#'
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param p Dimension of observation vecor.
#' @param g Number of multivariate normal classes.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @param paralist A list containing the required parameters \eqn{(\pi, \mu, \Sigma)}.
#' @param clust An n-dimensional vector of class partition.
#' @return
#' \item{errval}{a value of error rate}
#' @export
#' @details
#' The  error rate of the Bayes rule for a g-class Gaussian mixture model is given by
#' \deqn{
#' err(y;\theta)=1-\sum_{i=1}^g\pi_i Pr\{R(y;\theta)=i\mid Z \in C_i\}.
#' }
#' Here, we write
#' \deqn{
#' Pr\{R(y;\theta) \in C_i\mid Z\in C_i\}=\frac{\sum_{j=1}^nI_{C_i}(z_j)Q[z_j,R(y;\theta) ]}{\sum_{j=1}^nI_{C_i}(z_j)},
#' }
#' where \eqn{Q[u,v]=1} if \eqn{u=v} and \eqn{Q[u,v]=0} otherwise, and \eqn{I_{C_i}(z_j)} is an indicator function for the \eqn{i}th class.
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
#' zm<-dat$clust
#' zm[m==1]<-NA
#' inits<-initialvalue(g=4,zm=zm,dat=dat$Y)
#' \donttest{
#' fit_pc<-gmmsslm(dat=dat$Y,zm=zm,pi=inits$pi,mu=inits$mu,sigma=inits$sigma,xi=xi,type='full')
#' parlist<-paraextract(fit_pc)
#' erate(dat=dat$Y,p=3,g=4,paralist=parlist,clust=dat$clust)
#' }
#'
erate<-function(dat,p,g,pi=NULL,mu =NULL,sigma=NULL,paralist=NULL,clust){
  if (!is.null(paralist)) {
    pi <- paralist$pi
    mu <- paralist$mu
    sigma <- paralist$sigma
  }else{
    paralist$pi<-pi
    paralist$mu<-mu
    paralist$sigma<-sigma
  }
  n <- dim(dat)[1]
  est_clust<-bayesclassifier(dat,p=p,g=g,paralist=paralist)
  prob<-sapply(1:g,function(j)(1-sum(est_clust[clust==j]==j)/sum(clust==j)))
  errval<-sum(paralist$pi*(prob))
  return(errval)
}
