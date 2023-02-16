#' Normalize log-probability
#'
#' Normalize log-probability.
#' @param x A variable vector.
#' @return
#' \item{val}{A normalize log probability of variable vector.}
#' @export

normalise_logprob <- function(x){
  x <- exp(x-max(x))
  val<-x/sum(x)
  return(val)
}
