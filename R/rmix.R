#' Normal mixture model generator.
#'
#' Generate random observations from the normal mixture distributions.
#' @param  n Number of observations.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @return
#' \item{Y}{An \eqn{n\times p} numeric matrix with samples drawn in rows.}
#' \item{Z}{ An \eqn{n\times g} numeric matrix; each row represents zero-one indicator variables defining the known class of origin of each.}
#' \item{clust}{An n-dimensional vector of class partition.}
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
#' @import mvtnorm
#' @export
rmix <- function(n,pi,mu,sigma){
g=length(pi)
ncov=ifelse(is.na(dim(sigma)[3]),1,2)
if(ncov==1){
  nn<-table(sample(1:g, n, replace = TRUE, prob = pi))
  p=dim(mu)[1]
  X=NULL
  for(j in 1:g){
    if (nn[j] > 0){
      X1<- mvtnorm::rmvnorm(nn[j], mean=mu[,j], sigma=as.matrix(sigma)) # class 1 :a random number generator for the multivariate normal distribution
      X=rbind(X,X1)
    } else if (nn[j]==0){
      X1=matrix(NA,nn[j],p)
      X=rbind(X,X1)
    }
  }  }else{
    nn<-table(sample(1:g, n, replace = TRUE, prob = pi))
    p=dim(mu)[1]
    X=NULL
    for(j in 1:g){
      if (nn[j] > 0){
        X1<- mvtnorm::rmvnorm(nn[j], mean=mu[,j], sigma=as.matrix(sigma[,,j])) # class 1 :a random number generator for the multivariate normal distribution
        X=rbind(X,X1)
      } else if (nn[j]==0){
        X1=matrix(NA,nn[j],p)
        X=rbind(X,X1)
      }
    }
  }
y <- rep(1:g,nn)
rperm <- sample(n) #re-arrange order
y <- y[rperm]
X <- X[rperm,,drop=FALSE]
clust<-y
Y<-X
Z<-makelabelmatrix(y)
dat <- list(Y=Y,Z=Z,clust=clust)

return(dat)
}
