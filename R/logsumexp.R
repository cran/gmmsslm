#' log summation of exponential function
#'
#' log summation of exponential variable vector.
#' @param x A variable vector.
#' @return
#' \item{val}{log summation of exponential variable vector.}
#' @export
#'
logsumexp <- function(x){
  log(sum(exp(x - max(x)))) + max(x)
}
