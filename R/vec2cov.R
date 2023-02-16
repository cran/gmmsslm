#' Transform a vector into a matrix
#'
#' Transform a vector into a matrix i.e., Sigma=R^T*R
#' @param par A vector representing a variance matrix
#' @return sigma A variance matrix
#' @export
#' @details The variance matrix is decomposed by computing the Choleski factorization of a real symmetric positive-definite square matrix.
#' Then, storing the upper triangular factor of the Choleski decomposition into a vector.



vec2cov <- function(par){
  q <- length(par)
  p <- (-1+sqrt(1+4*1*2*q))/2
  R <- matrix(0, p,p)
  diag_elements <- par[1:p]
  upper_elements <-par[-(1:p)]
  diag(R) <- exp(diag_elements)
  if (any(is.infinite(R) | is.nan(R))) stop('Variances infinite or NaN')
  R[upper.tri(R)] <- upper_elements
  sigma <- t(R) %*% R
  if(any(eigen(sigma)$values <= 0)) stop('Conversion failed (negative eigenvalues)')
  return(sigma)
}
