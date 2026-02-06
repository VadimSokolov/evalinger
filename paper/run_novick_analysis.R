#!/usr/bin/env Rscript
# =============================================================================
# run_novick_analysis.R
#
# Real-data demonstration using Novick & Grizzle (1965) ulcer surgery trial.
#
# The data come from a 19-hospital VA cooperative study comparing four
# operations for duodenal ulcer (A, B, C, D), each with 100 patients.
# Outcomes were categorized as Excellent/Good, Fair/Poor, or Death.
#
# Death counts:  A = 7,  B = 1,  C = 1,  D = 3
#
# We treat the data as arriving sequentially (as in the original trial)
# and apply e-value monitoring to detect mortality differences.
#
# Outputs (saved to sims/ directory):
#   sims/novick_summary.rds          - summary table
#   sims/novick_eprocess.pdf         - e-process trajectories figure
#   sims/novick_confseq.pdf          - confidence sequence figure
#   sims/novick_platform.rds         - platform monitoring results
# =============================================================================

cat("=== Novick (1965) real-data analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Load evalinger package ---
# This script lives in evalinger/paper/; the package root is one level up.
pkg_dir <- if (file.exists(file.path("..", "DESCRIPTION"))) ".." else "."
cat("Loading evalinger from:", normalizePath(pkg_dir), "\n")
devtools::load_all(pkg_dir, quiet = TRUE)

sims_dir <- "sims"
if (!dir.exists(sims_dir)) dir.create(sims_dir, recursive = TRUE)

# =============================================================================
# 1.  DATA
# =============================================================================
# Novick & Grizzle (1965), Table 2
# Binary outcome: Death (1) vs Survival (0)
set.seed(1965)  # Year of the paper

n_per_arm <- 100L
deaths <- c(A = 7L, B = 1L, C = 1L, D = 3L)

make_arm <- function(n_deaths, n_total) {
  sample(c(rep(1L, n_deaths), rep(0L, n_total - n_deaths)))
}

x_A <- make_arm(deaths["A"], n_per_arm)
x_B <- make_arm(deaths["B"], n_per_arm)
x_C <- make_arm(deaths["C"], n_per_arm)
x_D <- make_arm(deaths["D"], n_per_arm)

cat("Death counts:  A =", deaths["A"], " B =", deaths["B"],
    " C =", deaths["C"], " D =", deaths["D"], "\n")

# =============================================================================
# 2.  PAIRWISE E-PROCESS MONITORING
# =============================================================================
# Test whether each arm has a HIGHER death rate than Treatment B (best performer)
# H_0: p_k = p_B   vs   H_1: p_k > p_B
#
# Design alternative: 5 percentage-point difference in death rate
# (consistent with Novick's prior expectation of 2-3% baseline mortality)
alpha <- 0.025
design_p_T <- 0.06   # hypothesized inferior-arm death rate
design_p_C <- 0.01   # reference-arm death rate
lambda_star <- grow_lambda(design_p_T, design_p_C)
cat("GROW-optimal lambda (design: 6% vs 1%):", round(lambda_star, 3), "\n")
cat("Expected growth rate:", round(expected_growth_rate(lambda_star,
    design_p_T, design_p_C), 4), "\n")
cat("Expected stopping time (approx):",
    round(expected_stopping_time(lambda_star, design_p_T, design_p_C,
                                 alpha = alpha)), "patients per arm\n\n")

# Pairwise e-processes (testing each arm vs B for higher death rate)
ep_AB <- eprocess_binary(x_A, x_B, lambda = lambda_star, alpha = alpha)
ep_AC <- eprocess_binary(x_A, x_C, lambda = lambda_star, alpha = alpha)
ep_AD <- eprocess_binary(x_A, x_D, lambda = lambda_star, alpha = alpha)
ep_DB <- eprocess_binary(x_D, x_B, lambda = lambda_star, alpha = alpha)
ep_DC <- eprocess_binary(x_D, x_C, lambda = lambda_star, alpha = alpha)

# C vs B (no true difference -- both 1% death rate)
ep_CB <- eprocess_binary(x_C, x_B, lambda = lambda_star, alpha = alpha)

cat("--- Pairwise e-process results (n = 100 per arm) ---\n")
results <- data.frame(
  Comparison = c("A vs B", "A vs C", "A vs D", "D vs B", "D vs C", "C vs B"),
  Deaths_T   = c(deaths["A"], deaths["A"], deaths["A"],
                 deaths["D"], deaths["D"], deaths["C"]),
  Deaths_C   = c(deaths["B"], deaths["C"], deaths["D"],
                 deaths["B"], deaths["C"], deaths["B"]),
  Final_E    = round(c(
    max(ep_AB$evalue), max(ep_AC$evalue), max(ep_AD$evalue),
    max(ep_DB$evalue), max(ep_DC$evalue), max(ep_CB$evalue)
  ), 2),
  Log_E      = round(c(
    max(ep_AB$log_evalue), max(ep_AC$log_evalue), max(ep_AD$log_evalue),
    max(ep_DB$log_evalue), max(ep_DC$log_evalue), max(ep_CB$log_evalue)
  ), 3),
  AV_pvalue  = round(pmin(1, 1 / c(
    max(ep_AB$evalue), max(ep_AC$evalue), max(ep_AD$evalue),
    max(ep_DB$evalue), max(ep_DC$evalue), max(ep_CB$evalue)
  )), 4),
  Rejected   = c(ep_AB$rejected, ep_AC$rejected, ep_AD$rejected,
                 ep_DB$rejected, ep_DC$rejected, ep_CB$rejected),
  Novick_Cred = c(0.980, 0.980, 0.887, 0.813, 0.813, NA),
  stringsAsFactors = FALSE
)
print(results)
cat("\nThreshold: E >= ", round(1/alpha, 1), " (log E >= ",
    round(log(1/alpha), 3), ")\n\n")

# =============================================================================
# 3.  CONFIDENCE SEQUENCE FOR A vs B
# =============================================================================
cs_AB <- confseq_binary(x_A, x_B, alpha = 0.05, t_min = 10)
cat("Confidence sequence for delta = p_A - p_B (death rate):\n")
cat("  Final point estimate:", round(cs_AB$delta_hat[n_per_arm], 3), "\n")
cat("  Final 95% CS: [", round(cs_AB$lower[n_per_arm], 3), ",",
    round(cs_AB$upper[n_per_arm], 3), "]\n\n")

# =============================================================================
# 4.  PLATFORM MONITORING WITH E-BH
# =============================================================================
# Treat B as shared control; test A, C, D
look_times <- c(25L, 50L, 75L, 100L)
pm <- platform_monitor(
  K = 3,
  look_times = look_times,
  x_C = x_B,
  x_T_list = list(x_A, x_C, x_D),
  lambda_vec = rep(lambda_star, 3),
  alpha = alpha,
  fdr_alpha = 0.05
)

cat("--- Platform monitoring (B = control; arms A, C, D) ---\n")
pm$summary$arm_label <- c("A", "C", "D")
print(pm$summary[, c("arm_label", "individual_rejected",
                      "individual_reject_time", "ebh_rejected_final")])

cat("\ne-BH results at each look:\n")
for (j in seq_along(look_times)) {
  cat("  Look", j, "(n =", look_times[j], "per arm):",
      pm$ebh_results[[j]]$n_rejected, "rejected\n")
}

# =============================================================================
# 5.  FIGURES
# =============================================================================

# --- Figure: E-process trajectories ---
pdf(file.path(sims_dir, "novick_eprocess.pdf"), width = 7, height = 5)
nn <- seq_len(n_per_arm)
threshold <- log(1 / alpha)

# Raw log-e-value trajectories (staircase pattern shows when evidence accumulates)
log_e_AB <- ep_AB$log_evalue
log_e_AD <- ep_AD$log_evalue
log_e_DB <- ep_DB$log_evalue
log_e_CB <- ep_CB$log_evalue

ylim <- range(c(log_e_AB, log_e_AD, log_e_DB, log_e_CB, threshold, 0))

par(mar = c(4.5, 4.5, 2, 1))
plot(nn, log_e_AB, type = "s", col = "#D55E00", lwd = 2,
     xlab = "Patient pairs", ylab = expression(log~E[n]),
     ylim = ylim,
     main = "E-process trajectories: Novick (1965) death rates")
lines(nn, log_e_AD, type = "s", col = "#0072B2", lwd = 2)
lines(nn, log_e_DB, type = "s", col = "#009E73", lwd = 2)
lines(nn, log_e_CB, type = "s", col = "gray50", lwd = 1.5, lty = 2)
abline(h = threshold, lty = 3, col = "black", lwd = 1.5)
abline(h = 0, col = "gray80")

legend("topleft",
       legend = c(
         paste0("A vs B  (7% vs 1%)"),
         paste0("A vs D  (7% vs 3%)"),
         paste0("D vs B  (3% vs 1%)"),
         paste0("C vs B  (1% vs 1%)"),
         expression(paste("Threshold  ", log(1/alpha)))
       ),
       col = c("#D55E00", "#0072B2", "#009E73", "gray50", "black"),
       lwd = c(2, 2, 2, 1.5, 1.5),
       lty = c(1, 1, 1, 2, 3),
       cex = 0.8, bty = "n")
dev.off()
cat("Saved:", file.path(sims_dir, "novick_eprocess.pdf"), "\n")

# --- Figure: Confidence sequence for A vs B ---
pdf(file.path(sims_dir, "novick_confseq.pdf"), width = 7, height = 4.5)
par(mar = c(4.5, 4.5, 2, 1))

# Start plotting from t_min = 10
idx <- 10:n_per_arm
plot(idx, cs_AB$delta_hat[idx], type = "l", lwd = 2,
     ylim = c(min(cs_AB$lower[idx]), max(cs_AB$upper[idx])),
     xlab = "Patient pairs", ylab = expression(hat(delta) == hat(p)[A] - hat(p)[B]),
     main = "Confidence sequence: Treatment A vs B mortality difference")
polygon(c(idx, rev(idx)),
        c(cs_AB$lower[idx], rev(cs_AB$upper[idx])),
        col = adjustcolor("steelblue", 0.25), border = NA)
lines(idx, cs_AB$lower[idx], col = "steelblue", lty = 2)
lines(idx, cs_AB$upper[idx], col = "steelblue", lty = 2)
lines(idx, cs_AB$delta_hat[idx], lwd = 2)
abline(h = 0, lty = 3, col = "gray50")
legend("topright",
       legend = c("Point estimate", "95% confidence sequence", "Null (no difference)"),
       col = c("black", "steelblue", "gray50"),
       lwd = c(2, 1, 1), lty = c(1, 2, 3), cex = 0.8, bty = "n")
dev.off()
cat("Saved:", file.path(sims_dir, "novick_confseq.pdf"), "\n")

# =============================================================================
# 6.  SAVE RESULTS
# =============================================================================
g_star <- expected_growth_rate(lambda_star, design_p_T, design_p_C)
est_stop <- expected_stopping_time(lambda_star, design_p_T, design_p_C,
                                    alpha = alpha)

novick_summary <- list(
  deaths = deaths,
  n_per_arm = n_per_arm,
  alpha = alpha,
  lambda = lambda_star,
  growth_rate = g_star,
  expected_stopping_time = est_stop,
  design_alternative = c(p_T = design_p_T, p_C = design_p_C),
  results = results,
  cs_final = list(
    estimate = cs_AB$delta_hat[n_per_arm],
    lower = cs_AB$lower[n_per_arm],
    upper = cs_AB$upper[n_per_arm]
  ),
  platform = pm$summary
)

saveRDS(novick_summary, file.path(sims_dir, "novick_summary.rds"))
saveRDS(pm, file.path(sims_dir, "novick_platform.rds"))

cat("\nAll Novick analysis outputs saved to", sims_dir, "\n")
cat("End time:", format(Sys.time()), "\n")
