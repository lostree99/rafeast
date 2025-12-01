#' Randomized-Accelerated FEAST
#' @export
rafeast_sparse <- function(A, lambda_min, lambda_max, m0, Nc=8, max_iter=10) {
  
  # --- Phase 1: Randomized Warmstart ---
  n <- nrow(A)
  Omega <- matrix(rnorm(n * (m0 + 10)), n, m0) 
  Y <- A %*% Omega
  qr_res <- qr(Y)
  Q_warm <- qr.Q(qr_res)
  Q_warm <- as.matrix(Q_warm[, 1:m0])

  # --- Phase 2: Refinement ---
  res <- feast_sparse(A, m0, lambda_min, lambda_max, 
                      q0 = Q_warm, 
                      Nc = Nc,            # <--- Pass this
                      max_iter = max_iter) # <--- Pass this
  
  return(res)
}