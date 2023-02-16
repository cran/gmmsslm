#' Transform a variance matrix into a vector
#'
#' Transform a variance matrix into a vector i.e., Sigma=R^T*R
#' @param sigma A \eqn{p\times p} variance matrix
#' @return par A vector representing a variance matrix
#' @export
#' @details The variance matrix is decomposed by computing the Choleski factorization of a real symmetric positive-definite square matrix.
#' Then, storing the upper triangular factor of the Choleski decomposition into a vector.

cov2vec <- function(sigma){
  R <- chol(sigma)
  upper_elements <- R[upper.tri(R, diag=FALSE)]
  diag_elements <- diag(R)
  par<-c(log(diag_elements), upper_elements)
  return(par)
}
