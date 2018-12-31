#### conjugate.R ####

#' Conjugate Normal - Inverse Gamma Update
#'
#' Generates a value from the posterior distribution in the case where there
#' is a multivariate normal likelihood and an inverse gamma prior.\cr \cr
#' Argument model:\cr
#' y ~ Nn(mu, sig2*R)\cr
#' sig2 ~ IG(a,b)\cr
#' @param y vector of values at the likelihood level.
#' @param a prior shape value for inverse gamma.
#' @param b prior SCALE value for inverse gamma.
#' @param mu mean vector for multivariate normal.
#' @param R correlation matrix for multivariate normal.
#' @param R_inv scaled precision matrix, alternative specification to R.
#' @export

update_normal_invgamma <- function(y, a, b, mu, R, R_inv) {
  if(missing(R)){
    qf <- Mahalanobis(y, mu, prec = R_inv)
  } else if(missing(R_inv)) {
    qf <- Mahalanobis(y, mu, cov = R)
  } else {
    stop("Provide either R or R_inv, but not both")
  }
  return(1 / rgamma(1, .5 * length(y) + a, rate = .5 * qf + b))
}


#' Conjugate Multivariate Normal - Multivariate Normal Update
#'
#' Generates a value from the posterior distribution in the case where
#' there is a multivariate normal likelihood and a multivariate normal prior.\cr \cr
#' Argument Model:\cr
#' y ~ Nn(X*beta, Sig)\cr
#' beta ~ Np(mu, V)\cr
#' @param y vector of values at the likelihood level.
#' @param X fixed design matrix in likelihood.
#' @param mu prior mean vector.
#' @param Sig,Sig_inv likelihood covariance/precision matrix.
#' @param V,V_inv prior covariance/precision matrix.
#' @export

update_normal_normal <- function(y, X, mu, Sig, V, Sig_inv, V_inv) {
  z <- rnorm(ncol(X))
  if(missing(Sig_inv)){
    if(missing(V_inv)){
      return(as.numeric(.update_nn_C(z, y, X, mu, Sig, FALSE, V, FALSE)))
    } else if(missing(V)){
      return(as.numeric(.update_nn_C(z, y, X, mu, Sig, FALSE, V_inv, TRUE)))
    } else {
      stop("Provide either V or V_inv, but not both")
    }
  } else if(missing(Sig)){
    if(missing(V_inv)){
      return(as.numeric(.update_nn_C(z, y, X, mu, Sig_inv, TRUE, V, FALSE)))
    } else if(missing(V)){
      return(as.numeric(.update_nn_C(z, y, X, mu, Sig_inv, TRUE, V_inv, TRUE)))
    } else {
      stop("Provide either V or V_inv, but not both")
    }
  } else {
    stop("Provide either Sig or Sig_inv, but not both")
  }
}

#' Gaussian Process Update
#'
#' Returns the updated parameters for a gaussian process 
#' given evaluation points, observed data, a mean function, and a covariance function.\cr \cr
#' Argument model:\cr
#' y(time) ~ GP(mnfun(x),covfun(x1-x2))
#'
#' @param x,y coordinate vectors for the observed data points.
#' @param time locations to evaluate the curve at.
#' @param mnfun mean function that takes a single argument.
#'   If using a constant like 0 use function(z)\{0\}.
#' @param covfun covariance function that takes a single argument that represents
#'   an absolute distance. Distances are internally calculated using fields::rdist.
#' @return a vector corresponding to the time variable that represents the updated
#'   mean vector.
#' @export

update_gaussian_process <- function(x, y, time, mnfun, covfun, get_var = FALSE) {
  R <- covfun(fields::rdist(c(time,x)))
  m <- length(time)
  f <- length(x)
  
  R11 <- R[1:m,1:m]
  R22 <- R[(m+1):(m+f),(m+1):(m+f)]
  R12 <- R[1:m,(m+1):(m+f)]
  
  mu1 <- mnfun(time)
  mu2 <- mnfun(x)
  up_mean <- .update_gp_mean_C(y, mu1, mu2, R12, R22)
  if(!get_var) return(as.numeric(up_mean))
  up_var <- .update_gp_var_C(R11, R12, R22)
  return(list(up_mean = as.numeric(up_mean), up_var = up_var))
}
