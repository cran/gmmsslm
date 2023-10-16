#' Plot Missingness Mechanism and Boxplot
#'
#' This function plots the smoothed values of `-log(entropy)` against the missingness mechanism
#' and a boxplot of entropy for labeled vs. unlabeled observations.
#'
#' @param dat An \eqn{n\times p} matrix where each row represents an individual observation
#' @param g Number of multivariate normal classes.
#' @param parlist A list containing the required parameters \eqn{(\pi, \mu, \Sigma)}.
#' @param zm An n-dimensional vector containing the class labels including the missing-label denoted as NA.
#' @param bandwidth Bandwidth for kernel smoothing. Default is 5.
#' @param range.x Range for x values. Default is c(0, 5).
#' @param kernel Kernel type for smoothing. Default is 'normal'.
#' @param ylim The y-axis limits in the form of c(ylim[1], ylim[2]). Default is NULL.
#'
#' @return A plot.
#' @export
#' @importFrom graphics boxplot par
#'
plot_missingness <- function(dat, g, parlist, zm, bandwidth = 5, range.x = c(0, 5), ylim=NULL ,kernel = 'normal') {
  # Get dimensions of dat
  n <- dim(dat)[1]
  p <- dim(dat)[2]
  zm<-ifelse(is.na(zm),1,0)

  # Get entropy using the provided function (assuming it's available in your environment)
  mlfitentropy <- get_entropy(dat = dat, n = n, p = p, g = g, mu = parlist$mu, sigma = parlist$sigma, pi = parlist$pi)

  # Set up a 1x2 plot layout
  par(mfrow=c(1,2))
  par(cex.main=1.25, cex.lab=1.25, cex.axis=1.25)

  # Boxplot
  boxplot(mlfitentropy ~ ifelse(zm==0, 'Labeled','Unlabeled'), main='(a) Boxplot', xlab='Observations', ylab='Entropy', outline=FALSE)

  # Compute -log(entropy)
  logentropy <- -log(mlfitentropy)

  # Kernel smoothing
  kfit <- ksmooth(logentropy, zm, range.x = range.x, kernel = kernel, bandwidth = bandwidth)

  # Plot
  plot(kfit$x, kfit$y, xlim = range.x, ylim =ylim, type = 'l',
       main = '(b) Estimated missingness mechanism', xlab = '-log(entropy)', ylab = 'Probability of missing label')
}
