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

// [[Rcpp::export(".dmnorm_C")]]
arma::vec dmnorm(arma::mat x,  arma::rowvec mu,  arma::mat covprec, bool is_prec, bool unnorm = false, bool log = false) {
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

// [[Rcpp::export(".rmnorm_C")]]
arma::vec rmnorm(arma::vec z, arma::vec mu, arma::mat covprec, bool is_prec) {
  if (is_prec) {
    return mu + solve(chol(covprec), z);
  } else {
    return mu + chol(covprec).t() * z;
  }
}

// [[Rcpp::export(".update_nn_C")]]
arma::vec update_nn(arma::vec z, arma::vec y, arma::mat X, arma::vec mu, arma::mat Sig, bool inv_Sig, arma::mat V, bool inv_V){
  if (!inv_V) {
    V = inv_sympd(V);
  }
  
  if (!inv_Sig) {
    //Sig = inv_sympd(Sig);
    arma::mat vv = X.t() * solve(Sig, X) + V; 
    return rmnorm(z, solve(vv, (X.t() * solve(Sig, y) + V * mu)), vv, true);
  } else {
    arma::mat vv = X.t() * Sig * X + V; 
    return rmnorm(z, solve(vv, (X.t() * Sig * y + V * mu)), vv, true);
  }
}

// [[Rcpp::export(".update_gp_mean_C")]]
arma::vec update_gp_mean(arma::vec y, arma::vec mu1, arma::vec mu2, arma::mat R12, arma::mat R22){
  return mu1 + R12 * solve(R22, y - mu2); //This seems slow, but it is fastest.
}

// [[Rcpp::export(".update_gp_var_C")]]
arma::mat update_gp_var(arma::mat R11, arma::mat R12, arma::mat R22){
  return R11 - R12 * solve(R22, R12.t());
}

arma::vec rep4k(arma::vec u){
  int n = u.n_elem;
  double tmp3;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 1.0){
      out(i) = 0.0;
    } else if(ui >= 0.0){
      tmp3 = 1.0 - ui;
      out(i) = tmp3*tmp3*tmp3;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec rep3k(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 2.0){
      out(i) = 0.0;
    } else if(ui > 1.0){
      out(i) = 2.0 - 3.0*ui + 3.0/2.0*ui*ui - 1.0/4.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 3.0*ui - 9.0/2.0*ui*ui + 7.0/4.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec rep2k(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 3.0){
      out(i) = 0.0;
    } else if(ui > 2.0){
      out(i) = 9.0/2.0 - 9.0/2.0*ui + 3.0/2.0*ui*ui - 1.0/6.0*ui*ui*ui;
    } else if(ui > 1.0){
      out(i) = -3.0/2.0 + 9.0/2.0*ui - 3.0*ui*ui + 7.0/12.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 3.0/2.0*ui*ui - 11.0/12.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

arma::vec evenk(arma::vec u){
  int n = u.n_elem;
  double ui;
  arma::vec out(n);
  for(int i=0; i<n; i++){
    ui = u(i);
    if(ui > 4.0){
      out(i) = 0.0;
    } else if(ui > 3.0){
      out(i) = -1.0/6.0 * (ui-4)*(ui-4)*(ui-4);
    } else if(ui > 2.0){
      out(i) = -22.0/3.0 + 10.0*ui - 4.0*ui*ui + 1.0/2.0*ui*ui*ui;
    } else if(ui > 1.0){
      out(i) = 2.0/3.0 - 2.0*ui + 2.0*ui*ui - 1.0/2.0*ui*ui*ui;
    } else if(ui > 0.0){
      out(i) = 1.0/6.0*ui*ui*ui;
    } else {
      out(i) = 0.0;
    }
  }
  return out;
}

// [[Rcpp::export(".bs_even_C")]]
arma::mat bs_even(arma::vec time, int nk){
  arma::vec u;
  double dnk = nk * 1.0;
  u = dnk / time.max() * (time - time.min());
  int nbasis = nk + 3;
  int ni = nk - 2;
  int neven_basis = ni - 1;
  arma::mat out(time.n_elem, nbasis);
  out.col(0) = rep4k(u);
  out.col(1) = rep3k(u);
  out.col(2) = rep2k(u);
  out.col(nbasis - 1) = rep4k(4 + neven_basis - 1 - u);
  out.col(nbasis - 2) = rep3k(4 + neven_basis - 1 - u);
  out.col(nbasis - 3) = rep2k(4 + neven_basis - 1 - u);
  for(int i=0; i < neven_basis; i++){
    out.col(3 + i) = evenk(u - i);
  }
  return out;
}

