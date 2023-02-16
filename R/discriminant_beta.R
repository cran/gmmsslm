#' Discriminant function
#'
#' Discriminant function in the particular case of g=2 classes with an equal-covariance matrix
#' @param pi A g-dimensional vector for the initial values of the mixing proportions.
#' @param  mu A \eqn{p \times g} matrix for the initial values of the location parameters.
#' @param sigma A \eqn{p\times p} covariance matrix.
#' @return
#' \item{beta0}{An intercept of discriminant function}
#' \item{beta}{A coefficient of discriminant function}
#' @export
#' @details
#' Discriminant function in the particular case of g=2 classes with an equal-covariance matrix can be expressed
#' \deqn{d(y_i,\beta)=\beta_0+\beta_1 y_i,}
#' where \eqn{\beta_0=\log\frac{\pi_1}{\pi_2}-\frac{1}{2}\frac{\mu_1^2-\mu_2^2}{\sigma^2}} and \eqn{\beta_1=\frac{\mu_1-\mu_2}{\sigma^2}}.

#calculate beta0 and beta1
discriminant_beta <- function(pi, mu, sigma){
  isigma <- solve(sigma)
  beta0 <- log(pi[1])-log(pi[2])-1/2*(t(mu[,1]) %*% isigma %*% mu[,1] -t(mu[,2]) %*% isigma %*% mu[,2])
  beta0 <- as.numeric(beta0)
  beta <-  isigma %*% (mu[,1]-mu[,2])
  return(list(beta0=beta0, beta=beta))
}
