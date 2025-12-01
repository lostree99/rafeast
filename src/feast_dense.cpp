#include <Rcpp.h>
using namespace Rcpp;

// Bring in FEAST's own C interface / macros
extern "C" {
  #include "feast.h"
  #include "feast_dense.h"
}

// [[Rcpp::export]]
Rcpp::List feast_interval_dense(Rcpp::NumericMatrix A,
                                double emin,
                                double emax,
                                int m0) {

  int n = A.nrow();
  if (A.ncol() != n) {
    Rcpp::stop("A must be square.");
  }

  // Copy A into column-major buffer for FEAST
  std::vector<double> a(n * n);
  std::copy(A.begin(), A.end(), a.begin());

  // FEAST parameters
  int fpm[128];
  feastinit(fpm);   // declared by feast.h

  double epsout = 0.0;
  int loop = 0;

  int m0_local = m0;
  std::vector<double> e(m0_local);
  std::vector<double> x(n * m0_local);
  std::vector<double> res(m0_local);
  int m = 0;
  int info = 0;

  char uplo = 'F';  // 'F' = full; 'U' if A stores only upper

  dfeast_syev(&uplo, &n,
              a.data(), &n,
              fpm, &epsout, &loop,
              &emin, &emax,
              &m0_local, e.data(), x.data(),
              &m, res.data(), &info);

  if (info != 0) {
    Rcpp::stop("dfeast_syev failed with info = %d", info);
  }

  // Trim to m eigenpairs
  NumericVector evals(m);
  for (int i = 0; i < m; ++i) evals[i] = e[i];

  NumericMatrix evecs(n, m);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      evecs(i, j) = x[i + j * n];
    }
  }

  NumericVector res_R(m);
  for (int i = 0; i < m; ++i) res_R[i] = res[i];

  return List::create(
    _["values"]     = evals,
    _["vectors"]    = evecs,
    _["m"]          = m,
    _["epsout"]     = epsout,
    _["loop"]       = loop,
    _["residuals"]  = res_R
  );
}
