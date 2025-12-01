#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
    void feastinit_(int *fpm);

    void dfeast_scsrev_(
        const char *uplo, const int *n, const double *a, const int *ia, const int *ja,
        int *fpm, double *epsout, int *loop, const double *emin, const double *emax,
        const int *m0, double *e, double *x, int *m, double *res, int *info
    );
}

// [[Rcpp::export]]
List feast_sparse(S4 A, int m0, double emin, double emax, 
                  Nullable<NumericMatrix> q0 = R_NilValue,
                  int Nc = 8, int max_iter = 20) {
    
    // --- 1. Extract Matrix ---
    IntegerVector p = A.slot("p");
    IntegerVector i = A.slot("i");
    NumericVector x = A.slot("x");
    IntegerVector dim = A.slot("Dim");
    
    int n = dim[0];
    int nnz = x.length();
    
    // --- 2. Prepare Indices ---
    std::vector<int> ia(n + 1);
    std::vector<int> ja(nnz);
    for(int k = 0; k <= n; k++) ia[k] = p[k] + 1;
    for(int k = 0; k < nnz; k++) ja[k] = i[k] + 1;
    
    // --- 3. Prepare Outputs ---
    NumericVector eval(m0);
    NumericMatrix evec(n, m0);
    NumericVector res(m0);
    int m_found = 0;
    
    // --- 4. Initialize Parameters ---
    int fpm[64];
    feastinit_(fpm); // Safe init
    
    fpm[0] = 0;         // Runtime report (0=no print)
    fpm[1] = Nc;        // <--- USE USER INPUT FOR CONTOUR POINTS
    fpm[3] = max_iter;  // <--- USE USER INPUT FOR ITERATIONS
    
    // --- 5. Handle Warmstart ---
    if (q0.isNotNull()) {
        NumericMatrix q_init(q0);
        if (q_init.nrow() != n || q_init.ncol() != m0) {
            stop("Warmstart dimension mismatch");
        }
        std::copy(q_init.begin(), q_init.end(), evec.begin());
        fpm[5] = 1; // Tell FEAST to use the provided subspace
    }
    
    // --- 6. Call FEAST ---
    char uplo = 'U'; // Correct setting for sparse
    double epsout = 0.0;
    int loop = 0;
    int info = 0;
    
    dfeast_scsrev_(
        &uplo, &n, x.begin(), ia.data(), ja.data(), fpm, &epsout, &loop, 
        &emin, &emax, &m0, eval.begin(), evec.begin(), &m_found, res.begin(), &info
    );
    
    if (info != 0 && info != 108) { 
        warning("FEAST returned info code: %d", info);
    }

    // --- 7. Return Results ---
    if (m_found < m0) {
        return List::create(
            Named("values") = head(eval, m_found),
            Named("vectors") = evec(Range(0, n-1), Range(0, m_found-1)),
            Named("nconv") = m_found
        );
    }
    
    return List::create(
        Named("values") = eval,
        Named("vectors") = evec,
        Named("nconv") = m_found
    );
}