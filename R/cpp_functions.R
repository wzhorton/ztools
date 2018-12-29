#### cpp_functions.R ####

#' Determinant of Symmetric Positive Definite Matrix
#' 
#' Computes the (log) determinant of a SymPD matrix using the Cholesky factorization.
#' @useDynLib ztools
#' @importFrom Rcpp sourceCpp
#' @exportPattern "^[[:alpha:]]+"
#' @param x a numeric matrix.
#' @param log logical; TRUE returns the log determinant.
#' @export

det_sympd <- function(x, log = FALSE){
  .det_sympd_C(as.matrix(x), log)
}

#' Fast Mahalanobis Distance Computation
#' 
#' Quickly evaluates Mahalanobis distance given either a covariance or precision matrix.
#' @param x a vector or matrix. Multiple observations should correspond to rows.
#' @param center the mean or centering vector
#' @param cov,prec covariance or precision matrices. Only provide one.
#' @export

Mahalanobis <- function(x, center, cov, prec){
  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  }
  if(missing(prec)) {
    return(.Mahalanobis_C(x, center, cov, FALSE))
  } else if(missing(cov)) {
    return(.Mahalanobis_C(x, center, prec, TRUE))
  }
  stop("Provide either Precision or Covariance, but not both")
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
    return(.dmnorm_C(x, mu, cov, FALSE, unnorm, log))
  } else if (missing(cov)) {
    out <- return(.dmnorm_C(x, mu, prec, TRUE, unnorm, log))
  } else {
    stop("Provide either cov or prec, but not both")
  }
}

