#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
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

// [[Rcpp::export]]
arma::vec dmnorm(arma::mat x,  arma::rowvec mu,  arma::mat covprec, bool is_prec, bool unnorm = false, bool log = false) {
  arma::vec distval = Mahalanobis(x,  mu, covprec, is_prec);
  double logdet;
  double log2pi;
  if (unnorm) {
    logdet = 0;
    log2pi = 0;
  } else {
    logdet = sum(arma::log(arma::eig_sym(covprec)));
    log2pi = std::log(2.0 * M_PI);
  }
  if (is_prec) {
    logdet = -1*logdet;
  }
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2 );
  
  if (log) {
    return(logretval);
  } else {
    return(exp(logretval));
  }
}
