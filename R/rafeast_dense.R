#' RA-FEAST on a dense SPD matrix
#'
#' @param A symmetric numeric matrix
#' @param lambda_min, lambda_max eigenvalue interval
#' @param m0 working subspace dimension (>= expected eigen count)
#' @param Nc number of contour points (default 16)
#' @param max_iter number of RA-FEAST iterations (default 8)
#' @return list with $values, $vectors, $m, $subspace_dim, $iter
#' @export
rafeast_dense <- function(A, lambda_min, lambda_max,
                          m0,
                          Nc = 16L,
                          max_iter = 8L) {
  if (!is.matrix(A)) A <- as.matrix(A)
  Nc <- as.integer(Nc)
  max_iter <- as.integer(max_iter)
  m0 <- as.integer(m0)

  rafeast_dense_cpp(A, lambda_min, lambda_max, m0, Nc, max_iter)
}
