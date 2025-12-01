#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# sweep_final_comprehensive.R (SMART RESUME VERSION)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(rafeast535)
  library(Matrix)
  library(igraph)
  library(RSpectra)
})

# --- 1. Job Index Setup ---
args <- commandArgs(trailingOnly = TRUE)
idx_arg <- if (length(args) >= 1) as.integer(args[1]) else NA_integer_
idx_env <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA)))
idx <- if (!is.na(idx_arg)) idx_arg else if (!is.na(idx_env)) idx_env else 1L

# --- 2. Check for Existing Result (The Fix) ---
dir.create("results_final", showWarnings = FALSE)
outfile <- sprintf("results_final/res_%05d.csv", idx)

# If file exists, check if it is valid (not empty, has data)
if (file.exists(outfile)) {
  tryCatch({
    existing_data <- read.csv(outfile)
    if (nrow(existing_data) > 0 && !any(is.na(existing_data$time_total))) {
      cat(sprintf("Job %d: Valid result exists. Skipping.\n", idx))
      quit(save = "no", status = 0)
    }
  }, error = function(e) {
    # If read fails, file is corrupt. We will overwrite it.
    cat(sprintf("Job %d: Found corrupt file. Re-running.\n", idx))
  })
}

# --- 3. The Master Grid (Same as before) ---
n_vals <- c(2000L, 4000L, 8000L, 16000L, 24000L, 32000L)
seeds  <- 1:40

# Method Configs
methods_df <- data.frame(
  type = "Standard",
  Nc = 8L,
  max_iter = 20L,
  label = "Standard_FEAST"
)

raf_grid <- expand.grid(
  Nc = c(2L, 4L, 8L),
  max_iter = c(2L, 4L)
)
raf_grid$type <- "RA-FEAST"
raf_grid$label <- paste0("RA-FEAST_Nc", raf_grid$Nc, "_Iter", raf_grid$max_iter)

all_methods <- rbind(methods_df, raf_grid[, names(methods_df)])

job_grid <- expand.grid(
  n_index = seq_along(n_vals),
  seed    = seeds,
  method_idx = seq_len(nrow(all_methods))
)

total_jobs <- nrow(job_grid)
if (idx < 1L || idx > total_jobs) stop(sprintf("Index %d out of range", idx))

conf <- job_grid[idx, ]
n    <- n_vals[conf$n_index]
seed <- conf$seed
method_conf <- all_methods[conf$method_idx, ]

# FORCE SERIAL EXECUTION
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

cat(sprintf("Job %04d/%d: N=%d, Seed=%d, Method=%s\n", idx, total_jobs, n, seed, method_conf$label))

# --- 4. Generate Matrix ---
set.seed(seed)
radius <- sqrt(log(n) / n) * 1.5
g <- sample_grg(n, radius)
A_adj <- as_adjacency_matrix(g, type = "both", sparse = TRUE)
L <- Diagonal(n, rowSums(A_adj)) - A_adj

m0 <- 40L; lmin <- 0.001; lmax <- 5.0

# --- 5. Run RSpectra (Truth) ---
t_rs_start <- proc.time()
rs_out <- tryCatch(
  eigs_sym(L, k=m0+5, which="SM", opts=list(retvec=TRUE)), 
  error=function(e) NULL
)
time_rs <- (proc.time() - t_rs_start)[["elapsed"]]
truth_vals <- if (!is.null(rs_out) && !any(is.na(rs_out$values))) {
  rs_out$values[rs_out$values >= lmin & rs_out$values <= lmax]
} else { NA }

# --- 6. Run Method ---
time_phase1 <- 0; time_phase2 <- 0; vals_est <- NA; resid_norm <- NA

gc()
t_total_start <- proc.time()

if (method_conf$type == "Standard") {
  res <- tryCatch(
    rafeast535::feast_sparse(L, m0, lmin, lmax, q0=NULL, Nc=method_conf$Nc, max_iter=method_conf$max_iter),
    error = function(e) NULL
  )
  if (!is.null(res)) {
    time_phase2 <- (proc.time() - t_total_start)[["elapsed"]]
    vals_est <- sort(res$values)
  }
} else {
  t_p1_start <- proc.time()
  Omega <- matrix(rnorm(n * (m0 + 10)), n, m0)
  Q_warm <- qr.Q(qr(L %*% Omega))
  Q_warm <- as.matrix(Q_warm[, 1:m0])
  time_phase1 <- (proc.time() - t_p1_start)[["elapsed"]]
  
  t_p2_start <- proc.time()
  res <- tryCatch(
    rafeast535::feast_sparse(L, m0, lmin, lmax, q0=Q_warm, Nc=method_conf$Nc, max_iter=method_conf$max_iter),
    error = function(e) NULL
  )
  time_phase2 <- (proc.time() - t_p2_start)[["elapsed"]]
  if (!is.null(res)) vals_est <- sort(res$values)
}

time_total <- (proc.time() - t_total_start)[["elapsed"]]

# --- 7. Metrics & Save ---
err_vs_truth <- NA
if (!all(is.na(vals_est)) && !all(is.na(truth_vals))) {
  err_vs_truth <- max(sapply(vals_est, function(x) min(abs(x - truth_vals))))
}

if (!is.null(res) && !is.null(res$vectors)) {
  mid_idx <- max(1, length(vals_est) %/% 2)
  lam <- vals_est[mid_idx]
  vec <- res$vectors[, mid_idx]
  resid_norm <- max(abs((L %*% vec) - lam * vec))
}

out_df <- data.frame(
  idx = idx, n = n, seed = seed,
  method_label = method_conf$label, method_type = method_conf$type,
  Nc = method_conf$Nc, max_iter = method_conf$max_iter,
  time_rs = time_rs, time_total = time_total,
  time_phase1 = time_phase1, time_phase2 = time_phase2,
  error_vs_truth = err_vs_truth, residual_norm = resid_norm
)

write.csv(out_df, outfile, row.names = FALSE)
cat("Saved", outfile, "\n")