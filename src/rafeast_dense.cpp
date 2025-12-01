#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

extern "C" {
  void zgetrf_(int* m, int* n,
               std::complex<double>* a, int* lda,
               int* ipiv, int* info);

  void zgetrs_(const char* trans, int* n, int* nrhs,
               const std::complex<double>* a, int* lda,
               const int* ipiv,
               std::complex<double>* b, int* ldb,
               int* info);
}

Rcpp::List rafeast_dense_cpp(Rcpp::NumericMatrix A_r,
                             double lambda_min,
                             double lambda_max,
                             int m0,
                             int Nc,
                             int max_iter) {
  int n = A_r.nrow();
  if (A_r.ncol() != n) {
    stop("A must be square.");
  }
  if (m0 <= 0 || m0 > n) {
    stop("m0 must be in 1..n");
  }
  if (Nc <= 0) {
    stop("Nc must be positive");
  }
  if (max_iter <= 0) {
    stop("max_iter must be positive");
  }

  // Wrap A as Armadillo matrix (no copy)
  mat A(A_r.begin(), n, n, false);

  // Complex version of A
  cx_mat A_cx = cx_mat(A, zeros<mat>(n, n));

  // Initial random subspace Q (real), n x m0
  mat Q0 = randn<mat>(n, m0);
  mat Q, Rtmp;
  qr_econ(Q, Rtmp, Q0); // orthonormalize

  // Precompute circular contour nodes and weights
  typedef std::complex<double> cd;
  cd I(0.0, 1.0);

  double c    = 0.5 * (lambda_min + lambda_max);
  double half = 0.5 * (lambda_max - lambda_min);
  double r    = 0.55 * half;      // slightly wider than half-width

  std::vector<cd> z_nodes(Nc);
  std::vector<cd> w_nodes(Nc);

  for (int j = 0; j < Nc; ++j) {
    double theta = 2.0 * M_PI * j / (double)Nc;
    cd eitheta   = std::exp(I * theta);
    cd z         = cd(c, 0.0) + r * eitheta;
    z_nodes[j]   = z;

    // dz/dθ = i r e^{iθ}
    cd dz_dtheta = I * r * eitheta;
    // weight ≈ (1/Nc) * (dz/dθ) / (2π i)
    w_nodes[j] = dz_dtheta / (double)Nc / (I * 2.0 * M_PI);
  }

  // -------------------------------------------------------------
  // Precompute LU factorizations for each (z_j I - A)
  // -------------------------------------------------------------
  struct LUData {
    cx_mat LU;           // in-place LU of (zI - A)
    std::vector<int> ipiv;
  };

  std::vector<LUData> lu_factors(Nc);

  for (int j = 0; j < Nc; ++j) {
    cd z = z_nodes[j];

    // shifted = z I - A
    cx_mat shifted = -A_cx;
    for (int i = 0; i < n; ++i) {
      shifted(i, i) += z;
    }

    LUData lud;
    lud.LU   = shifted;           // copy
    lud.ipiv = std::vector<int>(n);

    int m   = n;
    int nn  = n;
    int lda = n;
    int info = 0;

    zgetrf_(&m, &nn,
            reinterpret_cast<std::complex<double>*>(lud.LU.memptr()),
            &lda,
            lud.ipiv.data(),
            &info);

    if (info != 0) {
      stop("zgetrf_ failed for contour node j=%d with info=%d", j, info);
    }

    lu_factors[j] = std::move(lud);
  }

  int it_used = 0;

  // -------------------------------------------------------------
  // RA-FEAST iterations using precomputed LU factors
  // -------------------------------------------------------------
  for (int it = 0; it < max_iter; ++it) {
    it_used = it + 1;

    // Lift Q to complex
    cx_mat Qc(Q, zeros<mat>(n, m0));

    // Filtered subspace
    cx_mat Qe(n, m0, fill::zeros);

    for (int j = 0; j < Nc; ++j) {
      const LUData& lud = lu_factors[j];

      // Solve (z_j I - A) X = Qc using precomputed LU factors
      cx_mat X = Qc;  // copy RHS

      int nn    = n;
      int nrhs  = m0;
      int lda   = n;
      int ldb   = n;
      int info  = 0;
      char trans = 'N';

      zgetrs_(&trans,
              &nn,
              &nrhs,
              reinterpret_cast<const std::complex<double>*>(lud.LU.memptr()),
              &lda,
              lud.ipiv.data(),
              reinterpret_cast<std::complex<double>*>(X.memptr()),
              &ldb,
              &info);

      if (info != 0) {
        stop("zgetrs_ failed at iteration %d contour node %d info=%d",
             it, j, info);
      }

      Qe += w_nodes[j] * X;
    }

    // Take real part and orthonormalize
    mat Q_real = real(Qe);
    qr_econ(Q, Rtmp, Q_real);
  }

  // -------------------------------------------------------------
  // Rayleigh-Ritz in current subspace Q
  // -------------------------------------------------------------
  mat T = Q.t() * A * Q;       // m0 x m0 symmetric
  vec theta;
  mat Y;
  eig_sym(theta, Y, T);        // eigenvalues ascending

  // Sort just in case
  uvec order = sort_index(theta);
  theta = theta(order);
  Y     = Y.cols(order);

  // Approximate eigenvectors in original space
  mat V = Q * Y;               // n x m0

  // Select eigenvalues in [lambda_min, lambda_max]
  std::vector<double> vals;
  std::vector<int> idx_keep;

  for (int i = 0; i < m0; ++i) {
    double ev = theta[i];
    if (ev >= lambda_min && ev <= lambda_max) {
      vals.push_back(ev);
      idx_keep.push_back(i);
    }
  }

  int m = (int)vals.size();

  NumericVector out_vals(m);
  for (int i = 0; i < m; ++i) {
    out_vals[i] = vals[i];
  }

  NumericMatrix out_vecs(n, m);
  for (int j = 0; j < m; ++j) {
    int col_idx = idx_keep[j];
    for (int i = 0; i < n; ++i) {
      out_vecs(i, j) = V(i, col_idx);
    }
  }

  return List::create(
    _["values"]       = out_vals,
    _["vectors"]      = out_vecs,
    _["m"]            = m,
    _["subspace_dim"] = m0,
    _["iter"]         = it_used
  );
}
