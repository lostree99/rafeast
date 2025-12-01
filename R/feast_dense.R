#' Dense FEAST on an SPD matrix within [emin, emax]
#'
#' @param A numeric matrix (dense, symmetric)
#' @param emin, emax eigenvalue interval
#' @param m0 working subspace dimension (>= expected eigen count)
#' @return list with $values, $vectors, $m, $residuals, etc.
#' @export
feast_band_dense <- function(A, emin, emax, m0) {
  if (!is.matrix(A)) A <- as.matrix(A)
  feast_interval_dense(A, emin, emax, m0)
}