#' Error rate of the Bayes rule for a g-class Gaussian mixture model
#'
#'Error rate of the Bayes rule for a g-class Gaussian mixture model
#'
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param parlist A list contains estimated \eqn{\pi}, \eqn{\mu} and \eqn{\sigma}.
#' @param clust An n-dimensional vector of class partition.
#' @return
#' \item{errval}{a value of error rate}
#' @details
#' The  error rate of the Bayes rule for a g-class Gaussian mixture model is given by
#' \deqn{
#' err(y;\theta)=1-\sum_{i=1}^g\pi_i Pr\{R(y;\theta)=i\mid Z \in C_i\}.
#' }
#' where
#' \deqn{
#' Pr\{R(y;\theta) \in C_i\mid Z\in C_i\}=\frac{\sum_{j=1}^nI_{C_i}(z_j)Q[z_j,R(y;\theta) ]}{\sum_{j=1}^nI_{C_i}(z_j)},
#' }
#' @export
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
#' erate(dat=dat$Y,parlist=fit_pc,clust=dat$clust)
#' }
#'
erate<-function(dat,parlist,clust){
  n<-dim(dat)[1]
  p<-dim(dat)[2]
  g<-length(parlist$parhat$pi)
  est_clust<-bayesclassifier(dat,n,p,g,parlist$parhat$pi,parlist$parhat$mu,parlist$parhat$sigma)
  prob<-sapply(1:g,function(j)(1-sum(est_clust[clust==j]==j)/sum(clust==j)))
  errval<-sum(pi*(prob))
  return(errval)
}
