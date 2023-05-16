#' @import methods
#' @export
# Define a new Rclass
#' @title gmmsslmFit Class
#' @description gmmsslmFit objects store the results of fitting Gaussian mixture models using the gmmsslm function.
#' @slot objective A numeric value representing the objective likelihood.
#' @slot ncov A numeric value representing the number of covariance matrices.
#' @slot convergence A numeric value representing the convergence value.
#' @slot iteration An integer value representing the number of iterations.
#' @slot obs A matrix containing the input data.
#' @slot m A logical vector representing label indicators.
#' @slot n An integer value representing the number of observations.
#' @slot p An integer value representing the number of variables.
#' @slot g An integer value representing the number of Gaussian components.
#' @slot type A character value representing the type of Gaussian mixture model.
#' @slot pi A numeric vector representing the mixing proportions.
#' @slot mu A matrix representing the location parameters.
#' @slot sigma An array representing the covariance matrix or list of covariance matrices.
#' @slot xi A numeric value representing the coefficient for a logistic function of the Shannon entropy.
#' @exportClass gmmsslmFit


#' @name gmmsslmFit-class
#' @title Class '"gmmsslmFit"'
#' @description An S4 class representing the result of fitting a Gaussian mixture model using gmmsslm()
#' @seealso gmmsslm
setClass("gmmsslmFit", slots = list(
  objective = "numeric",
  ncov = "numeric",
  convergence = "numeric",
  iteration = "integer",
  obs="matrix",
  m="logical",
  n = "integer",
  p = "integer",
  g = "integer",
  type = "character",
  pi="numeric",
  mu="matrix",
  sigma="array",
  xi="ANY"
))

#' @name gmmsslm
#' @title Fitting Gaussian mixture model to a complete classified dataset or an incomplete classified dataset with/without the missing-data mechanism.
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param zm An n-dimensional vector containing the class labels including the missing-label denoted as NA.
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix,or a list of g covariance matrices with dimension \eqn{p\times p \times g}.
#' It is assumed to fit the model with a common covariance matrix if \code{sigma} is a \eqn{p\times p} covariance matrix;
#' otherwise it is assumed to fit the model with unequal covariance matrices.
#' @param paralist A list containing the required parameters \eqn{(\pi, \mu, \Sigma)}.
#' @param xi A 2-dimensional vector containing the initial values of the coefficients in the logistic function of the Shannon entropy.
#' @param type Three types of Gaussian mixture models, 'ign' indicates fitting the model to a partially classified sample on the basis of the likelihood that ignores the missing label mechanism,
#' 'full' indicates fitting the model to a partially classified sample on the basis of the full likelihood, taking into account the missing-label mechanism,
#' and 'com' indicate fitting the model to a completed classified sample.
#' @param iter.max Maximum number of iterations allowed. Defaults to 500
#' @param eval.max Maximum number of evaluations of the objective function allowed. Defaults to 500
#' @param rel.tol Relative tolerance. Defaults to 1e-15
#' @param sing.tol Singular convergence tolerance; defaults to 1e-20.
#' @return A gmmsslmFit object containing the following slots:
#' \item{objective}{Value of objective likelihood}
#' \item{convergence}{Value of convergence}
#' \item{iteration}{Number of iterations}
#' \item{obs}{Input data matrix}
#' \item{n}{Number of observations}
#' \item{p}{Number of variables}
#' \item{g}{Number of Gaussian components}
#' \item{type}{Type of Gaussian mixture model}
#' \item{pi}{Estimated vector of the mixing proportions}
#' \item{mu}{Estimated matrix of the location parameters}
#' \item{sigma}{Estimated covariance matrix or list of covariance matrices}
#' \item{xi}{Estimated coefficient vector for a logistic function of the Shannon entropy}
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
#' fit_pc<-gmmsslm(dat=dat$Y,zm=zm,paralist=inits,xi=xi,type='full')
#' }
gmmsslm <-function(dat,zm,pi=NULL,mu=NULL,sigma=NULL,paralist=NULL,xi=NULL,type,iter.max=500,eval.max=500,rel.tol=1e-15,sing.tol=1e-15){
  if (!is.null(paralist)) {
    pi <- paralist$pi
    mu <- paralist$mu
    sigma <- paralist$sigma
  }
  Y<-dat
  g<-length(pi)
  p=dim(Y)[2]
  ncov=ifelse(is.na(dim(sigma)[3]),1,2)
  mm=is.na(zm)
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
  if(type=='full'){xihat<-parhat$xi}else{xihat<-xi}
  result <- new("gmmsslmFit",
                objective = fullopt$objective,
                ncov = ncov,
                convergence = fullopt$convergence,
                iteration = fullopt$iterations,
                obs=Y,
                m=mm,
                n=dim(Y)[1],
                p=p,
                g=g,
                type=type,
                pi=parhat$pi,
                mu=parhat$mu,
                sigma=parhat$sigma,
                xi=xihat)
  return(result)
}

#' @name summary
#' @title Summary method for gmmsslmFit objects
#' @description This function extracts summary information from a gmmsslmFit object, including objective value, ncov, convergence, iteration, and type.
#' @param object A gmmsslmFit object.
#' @rdname summary
#' @aliases summary,gmmsslmFit-method
#' @import methods
#' @method summary gmmsslmFit
#' @export
setGeneric("summary", function(object) {
  standardGeneric("summary")
})
setMethod("summary", "gmmsslmFit",
          function(object) {
            # Check the value of ncov and assign the corresponding description
            ncov_description <- if (object@ncov == 1) "A common covariance matrix" else "Unequal covariance matrices"

            # Table content
            table_list <- list(
              Likelihood = object@objective,
              VarianceStructure = ncov_description,
              Convergence = object@convergence,
              Iteration = object@iteration,
              TotalObservation = object@n,
              Dimension = object@p,
              ModelType = object@type
            )

            # Parameters content
            parameters_list <- list(
              pi = object@pi,
              mu = object@mu,
              sigma = object@sigma
            )

            # Add xi to the parameters_list if the type is 'full'
            if (object@type == "full") {
              parameters_list$xi <- object@xi
            }

            # Convert the table_list to a data frame and transpose it
            summary_table <- t(data.frame(table_list))

            # Print the summary table and parameters list
            cat("Table:\n")
            print(summary_table)
            cat("\nParameters:\n")
            print(parameters_list)
          })


#' @name paraextract
#' @title Extract parameter list from gmmsslmFit objects
#' @description This function extracts the parameters from a gmmsslmFit object, including p, g, pi, mu, and sigma.
#' @param object A gmmsslmFit object.
#' @rdname paraextract
#' @aliases paraextract,gmmsslmFit-method
#' @import methods
#' @method paraextract gmmsslmFit
#' @export
setGeneric("paraextract", function(object) {
  standardGeneric("paraextract")
})
setMethod("paraextract", "gmmsslmFit",
          function(object) {
            # Parameters content
            parameters_list <- list(
              pi = object@pi,
              mu = object@mu,
              sigma = object@sigma
            )

            # Add xi to the parameters_list if the type is 'full'
            if (object@type == "full") {
              parameters_list$xi <- object@xi
            }

          return(parameters_list)
          })

#' @name predict
#' @title Predict unclassified label
#' @description This function predicts unclassified label from a gmmsslmFit object.
#' @param object A gmmsslmFit object.
#' @rdname predict
#' @aliases predict,gmmsslmFit-method
#' @import methods
#' @method predict gmmsslmFit
#' @export
setGeneric("predict", function(object) {
  standardGeneric("predict")
})
setMethod("predict", "gmmsslmFit",
          function(object) {
            dat <- object@obs
            dat_ul <- dat[object@m,]
            prelist <- list(pi = object@pi, mu = object@mu, sigma = object@sigma)
            predlabel <- bayesclassifier(dat_ul, p = object@p, g = object@g, paralist = prelist)
            return(predlabel)
          })
