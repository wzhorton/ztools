#### distribution_functions.R ####


#' Inverse Gamma Density Function
#'
#' Computes the (log) density of the inverse gamma distribution using
#' either the scale or rate parametrization. The scale parametrization has
#' expected value of scale/(shape - 1)
#'
#' @param x vector of nonnegative values.
#' @param shape,scale shape and scale parameters. Must be stricly positive.
#' @param rate alternative way to specify scale.
#' @param log logical; if TRUE, density calculations are computed on the log scale.
#' @export

dinvgamma <- function(x, shape, rate, scale = 1 / rate, log = FALSE) {
  if (any(shape <= 0) || any(scale <= 0)) {
    stop("Shape and rate/scale must be positive")
  }
  if (!is.logical(log)) {
    stop("log needs to be logical")
  }
  
  a <- shape
  b <- scale
  out <- a * log(b) - lgamma(a) + (-a - 1) * log(x) - b / x
  out[is.nan(out)] <- -Inf
  
  if (log == FALSE) out <- exp(out)
  return(out)
}

#' Multivariate Normal Density Function
#'
#' Computes the (log) density of the multivariate normal distribution
#' using either the covariance or precision parametrization.
#'
#' @param x vector or matrix of values. Multiple observations should correspond to rows.
#' @param mu mean vector.
#' @param cov,prec covariance or precision matrix. Only specify one.
#' @param unnorm logical; TRUE ommits all terms not affected by x.
#' @param log logical; if TRUE, density calculations are computed on the log scale.
#' @export

dmnorm <- function(x, mu, cov, prec, unnorm = FALSE, log = FALSE) {
  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  }
  if (missing(prec)) {
    return(as.numeric(.dmnorm_C(x, mu, cov, FALSE, unnorm, log)))
  } else if (missing(cov)) {
    out <- return(as.numeric(.dmnorm_C(x, mu, prec, TRUE, unnorm, log)))
  } else {
    stop("Provide either cov or prec, but not both")
  }
}

#' Random Multivariate Normal Generator
#'
#' Generates a normally distributed vector given a mean vector and
#' either a covariance matrix or a precision matrix. Computationally, this method
#' has an advantage since an inverse is never taken.
#'
#' @param mu mean vector.
#' @param cov covariance matrix.
#' @param prec precision matrix.
#'
#' @return A randomly generated vector with the same length as mu.
#' @export

rmnorm <- function(mu, cov, prec) {
  if(missing(prec)) {
    return(as.numeric(.rmnorm_C(rnorm(length(mu)), mu, cov, FALSE)))
  } else if(missing(cov)) {
    return(as.numeric(.rmnorm_C(rnorm(length(mu)), mu, prec, TRUE)))
  }
  stop("Provide either Precision or Covariance, but not both")
}


