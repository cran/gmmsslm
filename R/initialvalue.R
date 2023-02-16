#' Initial values for ECM
#'
#' Inittial values for claculating the estimates based on solely on the classified features.
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param zm An n-dimensional vector containing the class labels including the missing-label denoted as NA.
#' @param g Number of multivariate normal classes.
#' @param ncov Options of structure of sigma matrix;  the default value is 2;
#'  \code{ncov} = 1 for a common covariance matrix;
#'  \code{ncov} = 2 for the unequal  covariance/scale matrices.
#' @return
#' \item{pi}{A g-dimensional  initial vector of the mixing proportions.}
#' \item{mu}{A initial  \eqn{p \times g} matrix of the location parameters.}
#' \item{sigma}{A \eqn{p\times p} covariance matrix if \code{ncov=1}, or a list of g covariance matrices with dimension \eqn{p\times p \times g} if \code{ncov=2}.}
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
#' inits<-initialvalue(g=4,zm=zm,dat=dat$Y,ncov=2)
#'
#'
#'
initialvalue <- function(dat,zm,g,ncov=2){
  Y<-dat
  #labelled indicators
  k<-zm[is.na(zm)==FALSE]
  Y<-as.matrix(Y)
  p<- ncol(Y)
  # labelled observations
  Y<-as.matrix(Y[is.na(zm)==FALSE,])
  # the number of labelled observation
  n <- length(k)
  nn=NULL
  pi=NULL
  mu <- matrix(0, p, g)
  sigmaa<-array(0,dim=c(p,p,g))
  for(i in 1:g){
    nn[i]<-sum(k==i)
    pi[i]<-nn[i]/n
    mu[,i]<-apply(as.matrix(Y[k==i,]),2,mean)
    mu[,i][is.nan(mu[,i])]=1
    sigmaa[,,i] <- t(Y[k==i,]-mu[,i])%*%((Y[k==i,])-mu[,i])
  }
  if(ncov==1){
    sigma<-apply(sigmaa,c(1,2),sum)
    sigma<-sigma/n
  }else{
    for(j in 1:g){
      sigmaa[,,j]<-sigmaa[,,j]/nn[j]
    }
    sigma<-sigmaa
  }
  return(list(pi=pi, mu=mu, sigma=sigma))
}
