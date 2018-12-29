#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".det_sympd_C")]] 
double det_sympd (arma::mat x, bool Log = false) {
  arma::mat cholx = chol(x);
  arma::vec y = log(cholx.diag());
  if( Log ) {
    return 2 * sum(y);
  } else {
    return exp(2 * sum(y));
  }
}

// [[Rcpp::export(".Mahalanobis_C")]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat covprec, bool is_prec){
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  if (is_prec) {
    return sum((x_cen * covprec) % x_cen, 1);
  } else {
    return sum((x_cen * inv_sympd(covprec)) % x_cen, 1);
  }
}

// [[Rcpp::export(.dmnorm_C)]]
arma::vec dmnorm(arma::vec x,  arma::rowvec mu,  arma::mat covprec, bool is_prec, bool unnorm = false, bool log = false) {
  arma::vec distval = Mahalanobis(x,  mu, covprec, is_prec);
  double logdet;
  double log2pi;
  if (unnorm) {
    logdet = 0;
    log2pi = 0;
  } else {
    logdet = det_sympd(covprec, true);
    log2pi = std::log(2.0 * M_PI);
  }
  if (is_prec) {
    logdet = -1*logdet;
  }
  arma::vec logretval = -((x.n_cols * log2pi + logdet + distval)/2);
  
  if (log) {
    return(logretval);
  } else {
    return(exp(logretval));
  }
}

// [[Rcpp::export(.rmnorm_C)]]
arma::vec rmnorm(arma::vec z, arma::vec mu, arma::mat covprec, bool is_prec) {
  if (is_prec) {
    return mu + solve(chol(covprec), z);
  } else {
    return mu + chol(covprec).t() * z;
  }
}
