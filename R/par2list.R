#' Transfer a vector into a list
#'
#' Transfer a vector into a list
#' @param par A vector with list information.
#' @param p Dimension of observation vecor.
#' @param g Number of multivariate normal classes.
#' @param ncov Options of structure of sigma matrix;  the default value is 2;
#'  \code{ncov} = 1 for a common covariance matrix that \code{sigma} is a \eqn{p\times p} matrix.
#'  \code{ncov} = 2 for the unequal  covariance/scale matrices that
#'  \code{sigma} represents a list of g matrices with dimension \eqn{p\times p \times g}.
#' @param type Three types to fit to the model, 'ign' indicates fitting the model on the basis of the likelihood that ignores the missing label mechanism,
#' 'full' indicates that the model to be fitted on the basis of the full likelihood, taking into account the missing-label mechanism,
#' and 'com' indicate that the model to be fitted to a completed classified sample.
#' @export
#' @return
#' \item{parlist}{Return a list including \code{mu}, \code{pi}, \code{sigma} and \code{xi}.}
#'
#'
#'
#'
par2list <- function(par, g, p,ncov=2,type=c('ign','full','com')){
  if(type=='com'){
    type='ign'
  }
  if(ncov==1){
    mu <- matrix(par[1:(p*2)], p, 2)
    par <- par[-(1:(2*p))]
    q <- p*(p+1)/2
    cholpars <- par[1:(q*1)]
    sigma <- matrix(0,p,p)
    sigma <- vec2cov(cholpars)
    par <- par[-(1:(q))]
    tpro <- par[1]
    pi <- vec2pro(tpro)
    parlist <- list(pi=pi, mu=mu, sigma=sigma)
    if (type=='full'){
      parlist$xi <- par[-1]
    }
      }else{
    mu <- matrix(par[1:(p*g)], p, g)
    par <- par[-(1:(p*g))]
    q <- p*(p+1)/2
    cholpars <- matrix(par[1:(q*g)], q, g)
    sigma <- array(0, dim=c(p,p,g))
    for (h in 1:g){
      sigma[,,h] <- vec2cov(cholpars[,h])
    }
    par <- par[-(1:(q*g))]
    tpro <- par[1:(g-1)]
    pi=vec2pro(tpro)
    par <- par[-(1:(g-1))]
    parlist <- list( pi=pi,mu=mu, sigma=sigma)
    if (type=='full'){
      parlist$xi <- par
    }
  }
  return(parlist)
}
