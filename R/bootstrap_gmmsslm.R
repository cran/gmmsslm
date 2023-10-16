
#' @title Bootstrap Analysis for gmmsslm
#' @description This file provides functions to perform bootstrap analysis on the results of the gmmsslm function.
#' @import methods

# Assuming the necessary structures and classes from gmmsslm.R are loaded and available.

#' @title Perform Bootstrap Analysis on gmmsslm
#' @description This function performs non-parametric bootstrap to assess the variability of the gmmsslm function outputs.
#' @param dat A matrix where each row represents an individual observation.
#' @param zm A matrix or data frame of labels corresponding to dat.
#' @param pi A numeric vector representing the mixing proportions.
#' @param mu A matrix representing the location parameters.
#' @param sigma An array representing the covariance matrix or list of covariance matrices.
#' @param paralist A list of parameters.
#' @param xi A numeric value representing the coefficient for a logistic function of the Shannon entropy.
#' @param type A character value indicating the type of Gaussian mixture model.
#' @param iter.max An integer indicating the maximum number of iterations.
#' @param eval.max An integer indicating the maximum number of evaluations.
#' @param rel.tol A numeric value indicating the relative tolerance.
#' @param sing.tol A numeric value indicating the singularity tolerance.
#' @param B An integer indicating the number of bootstrap samples.
#' @return A list containing mean and sd of bootstrap samples for pi, mu, sigma, and xi.

bootstrap_gmmsslm <- function(dat, zm, pi, mu, sigma, paralist, xi, type, iter.max=500, eval.max=500, rel.tol=1e-15, sing.tol=1e-15, B = 2000) {

  pi_values <- vector("list", B)
  mu_values <- array(0, dim = c(B, dim(mu)[1], dim(mu)[2]))
  sigma_values <- array(0, dim = c(B, dim(sigma)[1], dim(sigma)[2], dim(sigma)[3]))
  xi_values <- vector("numeric", B)

  for (i in 1:B) {
    indices <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
    resample_data <- dat[indices, ]
    resample_zm <- zm[indices, ]

    boot_result <- gmmsslm(dat = resample_data, zm = resample_zm, pi = pi, mu = mu, sigma = sigma, paralist = paralist, xi = xi, type = type, iter.max = iter.max, eval.max = eval.max, rel.tol = rel.tol, sing.tol = sing.tol)

    pi_values[[i]] <- boot_result@pi
    mu_values[i,,] <- boot_result@mu
    sigma_values[i,,,] <- boot_result@sigma
    xi_values[i] <- boot_result@xi
  }

  # Compute mean and sd for the bootstrap samples
  results <- list(
    pi_mean = mean(unlist(pi_values)),
    pi_sd = sd(unlist(pi_values)),
    mu_mean = apply(mu_values, c(2, 3), mean),
    mu_sd = apply(mu_values, c(2, 3), sd),
    sigma_mean = apply(sigma_values, c(2, 3, 4), mean),
    sigma_sd = apply(sigma_values, c(2, 3, 4), sd),
    xi_mean = mean(xi_values),
    xi_sd = sd(xi_values)
  )

  return(results)
}

# Additional utility functions can be added below if necessary.
