#' Transfer an informative vector to a probability vector
#'
#' Transfer an informative vector to a probability vector
#' @param vec An informative vector
#' @return pro A probability vector
#' @export

vec2pro <- function(vec){
  g <- length(vec)+1
  pro <- numeric(g)
  z <- numeric(g)
  for (h in 1:(g-1)){
    z[h] <- 1/(1+exp(-vec[h]-log(1/(g-h))))
  }
  pro[1] <- z[1]
  if (g > 2){
    for (h in 2:(g-1)){
      pro[h] <- (1-sum(pro[1:(h-1)]))*z[h]
    }
  }
  pro[g] <- 1-sum(pro[-g])
  return(pro)
}
