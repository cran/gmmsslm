#' Full log-likelihood function
#'
#' Full log-likelihood function with both terms of ignoring and missing
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param zm An n-dimensional vector containing the class labels including the missing-label denoted as NA.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @param xi A 2-dimensional vector containing the initial values of the coefficients in the logistic function of the Shannon entropy.
#' @details
#' The full log-likelihood function can be expressed as
#' \deqn{
#' \log L_{PC}^{({full})}(\boldsymbol{\Psi})=\log L_{PC}^{({ig})}(\theta)+\log L_{PC}^{({miss})}(\theta,\boldsymbol{\xi}),}
#' where\eqn{\log L_{PC}^{({ig})}(\theta)}is the log likelihood function formed ignoring the missing in the label of the unclassified features,
#' and \eqn{\log L_{PC}^{({miss})}(\theta,\boldsymbol{\xi})} is the log likelihood function formed on the basis of the missing-label indicator.
#' @return
#' \item{lk}{Log-likelihood value}
#' @export



#log likelihood L_pc^(full)
loglk_full<-function(dat,zm,pi,mu,sigma,xi){
  Y=dat
  lk_pc_ig<-loglk_ig(zm=zm,dat=dat,pi=pi,mu=mu,sigma=sigma)
  lk_pc_miss<-loglk_miss(zm=zm,dat=dat,xi=xi,pi=pi,mu=mu,sigma=sigma)
  lk_pc_full<-lk_pc_ig+lk_pc_miss
  return(lk_pc_full)
}
