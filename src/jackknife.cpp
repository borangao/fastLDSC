// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List run_jackknife_engine(const arma::mat& X,
                          const arma::vec& y,
                          int n_blocks,
                          double overlap_prop,
                          double scale) {
  int n_snps   = X.n_rows;
  int n_params = X.n_cols;

  // 1. block partition
  arma::uvec limits = arma::linspace<arma::uvec>(0, n_snps, n_blocks + 1);

  // 2. Pre-compute block-wise XtX and Xty
  arma::field<arma::mat> block_XtX(n_blocks);
  arma::field<arma::vec> block_Xty(n_blocks);

  arma::mat total_XtX(n_params, n_params, fill::zeros);
  arma::vec total_Xty(n_params, fill::zeros);

  for (int i = 0; i < n_blocks; i++) {
    int start = limits[i];
    int end   = (int)limits[i + 1] - 1;
    if (start > end) continue;

    arma::mat X_block = X.rows(start, end);
    arma::vec y_block = y.rows(start, end);

    arma::mat xtx = X_block.t() * X_block;
    arma::vec xty = X_block.t() * y_block;

    block_XtX(i) = xtx;
    block_Xty(i) = xty;

    total_XtX += xtx;
    total_Xty += xty;
  }

  // 3. full-data estimate
  arma::vec reg_tot = solve(total_XtX, total_Xty);

  double slope_tot     = reg_tot(0);
  double intercept_tot = reg_tot(1);

  double rhog_tot = slope_tot * scale;
  double rhoe_tot = 0.0;

  if (std::abs(overlap_prop) > 1e-9) {
    double total_pheno_corr = intercept_tot / overlap_prop;
    rhoe_tot = total_pheno_corr - rhog_tot;
  }

  // 4. Jackknife loop
  arma::mat pseudo_values(n_blocks, n_params, fill::zeros);
  arma::vec pseudo_rhoe(n_blocks, fill::zeros);

  for (int i = 0; i < n_blocks; i++) {
    arma::mat XtX_del = total_XtX - block_XtX(i);
    arma::vec Xty_del = total_Xty - block_Xty(i);

    arma::vec reg_del = solve(XtX_del, Xty_del);

    arma::vec pseudo = (double)n_blocks * reg_tot - ((double)n_blocks - 1.0) * reg_del;
    pseudo_values.row(i) = pseudo.t();

    double slope_del     = reg_del(0);
    double intercept_del = reg_del(1);

    double rhog_del = slope_del * scale;
    double rhoe_del = 0.0;

    if (std::abs(overlap_prop) > 1e-9) {
      double total_pheno_corr_del = intercept_del / overlap_prop;
      rhoe_del = total_pheno_corr_del - rhog_del;
    }

    pseudo_rhoe(i) = (double)n_blocks * rhoe_tot - ((double)n_blocks - 1.0) * rhoe_del;
  }

  return List::create(
    Named("reg")          = reg_tot,
    Named("pseudo")       = pseudo_values,
    Named("rho_e")        = rhoe_tot,
    Named("pseudo_rho_e") = pseudo_rhoe
  );
}

