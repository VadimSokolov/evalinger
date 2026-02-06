#!/usr/bin/env Rscript
# =============================================================================
# run_simulations.R
#
# Standalone simulation script for the paper:
#   "E-values for Adaptive Clinical Trials"
#
# Produces all tables, figures, and saved objects used by the paper.
# Run this ONCE (or when you want to update results), then render the paper.
#
# Outputs (saved to sims/ directory):
#   sims/sim_results.rds       - main comparison object (ecomparison)
#   sims/sensitivity.rds       - sensitivity analysis data frame
#   sims/hybrid.rds            - hybrid monitoring table
#   sims/concordance.rds       - concordance statistics
#   sims/design_params.rds     - design parameters and GROW values
#   sims/eprocess-paths.pdf    - e-process sample paths figure
#   sims/growth-landscape.pdf  - growth rate landscape figure
#   sims/hybrid-monitoring.pdf - hybrid monitoring figure
# =============================================================================

cat("=== E-values paper: running simulations ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Load evalinger package ---
# This script lives in evalinger/paper/; the package root is one level up.
pkg_dir <- file.path(dirname(getwd()), "")
if (file.exists(file.path("..", "DESCRIPTION"))) {
  pkg_dir <- ".."
} else if (file.exists("DESCRIPTION")) {
  pkg_dir <- "."
}
cat("Loading evalinger from:", normalizePath(pkg_dir), "\n")
devtools::load_all(pkg_dir, quiet = TRUE)

# --- Create output directory ---
out_dir <- "sims"
if (!dir.exists(out_dir)) dir.create(out_dir)

# --- Design parameters ---
p_C     <- 0.30
p_T_alt <- 0.45
Nmax    <- 200
n_looks <- 20
alpha   <- 0.025
nrep    <- 50000  # fast enough with vectorized implementation

lambda_star <- grow_lambda(p_T_alt, p_C)
g_star      <- expected_growth_rate(lambda_star, p_T_alt, p_C)
E_N         <- expected_stopping_time(lambda_star, p_T_alt, p_C, alpha)

design_params <- list(
  p_C = p_C, p_T_alt = p_T_alt, Nmax = Nmax, n_looks = n_looks,
  alpha = alpha, nrep = nrep,
  lambda_star = lambda_star, g_star = g_star, E_N = E_N
)
saveRDS(design_params, file.path(out_dir, "design_params.rds"))
cat(sprintf("GROW-optimal lambda: %.4f\n", lambda_star))
cat(sprintf("Expected growth rate: %.5f nats/pair\n", g_star))
cat(sprintf("Expected stopping time: %.0f pairs\n", E_N))

# =============================================================================
# 1. Main five-method comparison
# =============================================================================
cat("\n--- Running main comparison (nrep =", nrep, ") ---\n")
t0 <- Sys.time()
cmp <- simulate_comparison(
  p_C = p_C, p_T_alt = p_T_alt, Nmax = Nmax,
  n_looks = n_looks, alpha = alpha,
  nrep = nrep, seed = 42
)
cat("  Elapsed:", format(Sys.time() - t0), "\n")
saveRDS(cmp, file.path(out_dir, "sim_results.rds"))
cat("  Main comparison results:\n")
print(cmp$results)

# =============================================================================
# 2. Concordance between GS and e-value
# =============================================================================
cat("\n--- Computing concordance ---\n")
concordance <- list(
  both_rej = mean(cmp$raw$gs_obf$alt_rej & cmp$raw$evalue$alt_rej),
  neither  = mean(!cmp$raw$gs_obf$alt_rej & !cmp$raw$evalue$alt_rej),
  gs_only  = mean(cmp$raw$gs_obf$alt_rej & !cmp$raw$evalue$alt_rej),
  e_only   = mean(!cmp$raw$gs_obf$alt_rej & cmp$raw$evalue$alt_rej)
)
saveRDS(concordance, file.path(out_dir, "concordance.rds"))
cat(sprintf("  Both reject: %.1f%%\n", 100 * concordance$both_rej))
cat(sprintf("  Neither: %.1f%%\n", 100 * concordance$neither))
cat(sprintf("  GS only: %.1f%%\n", 100 * concordance$gs_only))
cat(sprintf("  E-value only: %.1f%%\n", 100 * concordance$e_only))

# =============================================================================
# 3. Sensitivity analysis over lambda
# =============================================================================
cat("\n--- Running sensitivity analysis ---\n")
lambda_grid <- c(0.10, 0.20, round(lambda_star, 2), 0.40, 0.50)
t0 <- Sys.time()
sens_results <- do.call(rbind, lapply(lambda_grid, function(lam) {
  cat(sprintf("  lambda = %.2f ... ", lam))
  res <- simulate_comparison(
    p_C = p_C, p_T_alt = p_T_alt, Nmax = Nmax,
    n_looks = n_looks, alpha = alpha,
    methods = "evalue", nrep = nrep, seed = 42, lambda = lam
  )
  cat("done\n")
  cbind(lambda = lam, res$results[, c("null_rej", "alt_rej",
                                       "avg_n_null", "avg_n_alt")])
}))
cat("  Elapsed:", format(Sys.time() - t0), "\n")
saveRDS(sens_results, file.path(out_dir, "sensitivity.rds"))

# =============================================================================
# 4. Hybrid monitoring (single trial)
# =============================================================================
cat("\n--- Running hybrid monitoring ---\n")
set.seed(42)
x_T_hybrid <- rbinom(Nmax, 1, p_T_alt)
x_C_hybrid <- rbinom(Nmax, 1, p_C)
look_times <- seq(Nmax / n_looks, Nmax, length.out = n_looks)

hm <- hybrid_monitor(x_T_hybrid, x_C_hybrid, look_times = look_times,
                     lambda = lambda_star, alpha = alpha,
                     p_T_design = p_T_alt, p_C_design = p_C,
                     gs_c = cmp$design$gs_c)
saveRDS(hm, file.path(out_dir, "hybrid.rds"))

# =============================================================================
# 5. Figure: E-process sample paths
# =============================================================================
cat("\n--- Generating e-process paths figure ---\n")
pdf(file.path(out_dir, "eprocess-paths.pdf"), width = 9, height = 4)
set.seed(123)
n_fig_paths <- 12
threshold <- log(1 / alpha)

par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))

# Null paths
plot(NULL, xlim = c(1, Nmax), ylim = c(-15, threshold + 2),
     xlab = "Observations per arm", ylab = expression(log(E[n])),
     main = expression("Under " * H[0]))
abline(h = threshold, col = "red", lty = 2, lwd = 2)
cols_null <- colorRampPalette(c("steelblue", "gray60"))(n_fig_paths)
for (i in seq_len(n_fig_paths)) {
  xT <- rbinom(Nmax, 1, p_C)
  xC <- rbinom(Nmax, 1, p_C)
  ep <- eprocess_binary(xT, xC, lambda = lambda_star, alpha = alpha)
  lines(seq_len(Nmax), ep$log_evalue, col = cols_null[i], lwd = 0.8)
}
legend("bottomleft", legend = c(expression(log(1/alpha)), "E-process paths"),
       col = c("red", "steelblue"), lty = c(2, 1), lwd = c(2, 1), bty = "n",
       cex = 0.8)

# Alternative paths
plot(NULL, xlim = c(1, Nmax), ylim = c(-10, 20),
     xlab = "Observations per arm", ylab = expression(log(E[n])),
     main = expression("Under " * H[1]))
abline(h = threshold, col = "red", lty = 2, lwd = 2)
cols_alt <- colorRampPalette(c("darkorange", "firebrick"))(n_fig_paths)
for (i in seq_len(n_fig_paths)) {
  xT <- rbinom(Nmax, 1, p_T_alt)
  xC <- rbinom(Nmax, 1, p_C)
  ep <- eprocess_binary(xT, xC, lambda = lambda_star, alpha = alpha)
  lines(seq_len(Nmax), ep$log_evalue, col = cols_alt[i], lwd = 0.8)
}
legend("topleft", legend = c(expression(log(1/alpha)), "E-process paths"),
       col = c("red", "darkorange"), lty = c(2, 1), lwd = c(2, 1), bty = "n",
       cex = 0.8)
dev.off()
cat("  Saved eprocess-paths.pdf\n")

# =============================================================================
# 6. Figure: Growth rate landscape
# =============================================================================
cat("\n--- Generating growth rate landscape figure ---\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  grid_df <- grow_lambda_grid(
    p_C = p_C,
    delta_grid = c(0.10, 0.15, 0.20),
    lambda_grid = seq(0.02, 0.70, by = 0.02),
    alpha = alpha
  )

  opt_df <- data.frame(
    delta = c(0.10, 0.15, 0.20),
    lam_star = sapply(c(0.10, 0.15, 0.20),
                      function(d) grow_lambda(p_C + d, p_C))
  )

  p_growth <- ggplot(grid_df, aes(x = lambda, y = growth_rate,
                                  color = factor(delta),
                                  group = factor(delta))) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
    geom_vline(data = opt_df,
               aes(xintercept = lam_star, color = factor(delta)),
               linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(
      values = c("0.1" = "steelblue", "0.15" = "darkorange", "0.2" = "firebrick"),
      labels = c(expression(delta == 0.10),
                 expression(delta == 0.15),
                 expression(delta == 0.20))
    ) +
    labs(x = expression(lambda), y = expression(g(lambda) ~ "(nats/pair)"),
         color = "Design alt.") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")

  ggsave(file.path(out_dir, "growth-landscape.pdf"), p_growth,
         width = 6, height = 3.5)
  cat("  Saved growth-landscape.pdf\n")
} else {
  cat("  ggplot2 not available; skipping growth landscape figure\n")
}

# =============================================================================
# 7. Figure: Hybrid monitoring
# =============================================================================
cat("\n--- Generating hybrid monitoring figure ---\n")
ep_hybrid <- eprocess_binary(x_T_hybrid, x_C_hybrid,
                             lambda = lambda_star, alpha = alpha)

pdf(file.path(out_dir, "hybrid-monitoring.pdf"), width = 7, height = 5)
plot_hybrid(ep_hybrid, look_times = look_times, Nmax = Nmax,
            alpha = alpha, x_T = x_T_hybrid, x_C = x_C_hybrid,
            gs_c = cmp$design$gs_c)
dev.off()
cat("  Saved hybrid-monitoring.pdf\n")

# =============================================================================
cat("\n=== All simulations complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("Outputs saved to:", normalizePath(out_dir), "\n")
