#' Fitting Gaussian mixture models
#'
#'Fitting Gaussian mixture model to a complete classified dataset or a incomplete classified dataset with/without the missing-data mechanism.
#'
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param zm An n-dimensional vector containing the class labels including the missing-label denoted as NA.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @param xi A 2-dimensional vector containing the initial values of the coefficients in the logistic function of the Shannon entropy.
#' @param type Three types of Gaussian mixture models, 'ign' indicates fitting the model to a partially classified sample on the basis of the likelihood that ignores the missing label mechanism,
#' 'full' indicates fitting the model to a partially classified sample on the basis of the full likelihood, taking into account the missing-label mechanism,
#' and 'com' indicate fitting the model to a completed classified sample.
#' @param iter.max Maximum number of iterations allowed. Defaults to 500
#' @param eval.max Maximum number of evaluations of the objective function allowed. Defaults to 500
#' @param rel.tol Relative tolerance. Defaults to 1e-15
#' @param sing.tol Singular convergence tolerance; defaults to 1e-20.
#' @return
#' \item{objective}{Value of objective likelihood}
#' \item{convergence}{Value of convergence}
#' \item{iteration}{Number of iteration}
#' \item{pi}{Estimated vector of the mixing proportions.}
#' \item{mu}{Estimated matrix of the location parameters.}
#' \item{sigma}{Estimated covariance matrix}
#' \item{xi}{Estimated  coefficient vector for a logistic function of the Shannon entropy}
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
#' }
#'

gmmsslm <-function(dat,zm,pi,mu,sigma,xi=NULL,type,iter.max=500,eval.max=500,rel.tol=1e-6,sing.tol=1e-20){
  Y<-dat
  g<-length(pi)
  p=dim(Y)[2]
  ncov=ifelse(is.na(dim(sigma)[3]),1,2)
  if(type=='com'){
    type='ign'
    if(any(is.na(zm))){
    stop('Missing labels exist in the completed classified sample')
    }
  }
  par<-list2par(pi=pi,mu=mu,sigma=sigma,xi=xi,p=p,g=g,type=type)
  fullopt<-nlminb(start=par,objective=neg_objective_function,gradient = NULL, hessian = NULL,g=g,zm=zm,dat=dat,type=type,control=list(iter.max=iter.max, eval.max=eval.max, rel.tol=rel.tol, sing.tol=sing.tol))
  #print(fullopt)
  parhat<-par2list(fullopt$par,g,p,type=type,ncov=ncov)
  parlist<-list(objective=fullopt$objective,ncov=ncov,convergence=fullopt$convergence,iteration=fullopt$iterations,parhat=parhat)
  return(parlist)
}
