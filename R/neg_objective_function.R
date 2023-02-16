#' Negative objective function for gmmssl
#'
#' Negative objective function for gmmssl
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param zm An n-dimensional vector of group partition including the missing-label, denoted as NA.
#' @param g Number of multivariate Gaussian groups.
#' @param par An informative vector including \code{mu}, \code{pi},\code{sigma} and \code{xi}.
#' @param ncov Options of structure of sigma matrix;  the default value is 2;
#'  \code{ncov} = 1 for a common covariance matrix;
#'  \code{ncov} = 2 for the unequal  covariance/scale matrices.
#' @param type Three types to fit to the model, 'ign' indicates fitting the model on the basis of the likelihood that ignores the missing label mechanism,
#' 'full' indicates that the model to be fitted on the basis of the full likelihood, taking into account the missing-label mechanism,
#' and 'com' indicate that the model to be fitted to a completed classified sample.
#'
#' @export
#'
#' @return
#'  \item{val}{Value of negatvie objective function.}
#'
neg_objective_function<-function(dat,zm,g,par,ncov=2,type=c('ign','full','com')){
  Y<-dat
  if(type=='com'){
    type='ign'
    if(any(is.na(zm))){
      stop('Missing labels exist in the completed classified sample')
    }
  }
  if(ncov==1){
    p <- ncol(Y)
    parlist <- par2list(par=par,g=g, p=p, type = type,ncov=1)
    pi <- parlist$pi
    mu <- parlist$mu
    sigma <- parlist$sigma
    if (type=='ign'){
      val <- loglk_ig(zm=zm, dat=dat, pi=pi, mu=mu, sigma=sigma)
    } else{
      xi=parlist$xi
      val=loglk_full(zm=zm, dat=dat,pi=pi,mu=mu,sigma=sigma,xi=xi)
    }

      }else{
    p <- ncol(Y)
    parlist <- par2list(par=par,g=g,p=p, type = type,ncov=ncov)
    pi=parlist$pi
    mu=parlist$mu
    sigma=parlist$sigma
    if (type=='ign'){
      val=loglk_ig(zm=zm, dat=dat, pi=pi, mu=mu, sigma=sigma)
    } else{
      xi=parlist$xi
      val=loglk_full(zm=zm,dat=dat,pi=pi,mu=mu,sigma=sigma,xi=xi)
    }
  }
  return(-val)
}
