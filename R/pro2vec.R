#' Transfer a probability vector into a vector
#'
#' Transfer a probability vector into an informative vector
#'
#' @param pro An propability vector
#' @return y An informative vector
#' @export

pro2vec <- function(pro){
  g <- length(pro)
  z <- numeric(g)
  z[1] <- pro[1]
  for (h in 2:g){
    z[h] <- pro[h]/(1-sum(pro[1:(h-1)]))
  }
  z <- z[-g]
  y <- log(z/(1-z))-log(1/(g-1:(g-1)))
  return(y)
}
