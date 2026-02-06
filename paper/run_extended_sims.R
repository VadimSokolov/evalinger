#!/usr/bin/env Rscript
# =============================================================================
# run_extended_sims.R
#
# Extended simulations for the paper:
#   "E-values for Adaptive Clinical Trials"
#
# Four new analyses:
#   1. Irregular/random look schedule: GS vs e-value
#   2. Parameter grid: vary p_C, delta, n_looks
#   3. Futility monitoring demonstration
#   4. RECOVERY-like large-trial simulation
#
# Outputs (saved to sims/ directory):
#   sims/irregular_comparison.rds  - irregular vs fixed look comparison
#   sims/parameter_grid.rds        - parameter grid results
#   sims/futility_demo.rds         - futility monitoring results
#   sims/recovery_sim.rds          - RECOVERY-like simulation results
#   sims/irregular-comparison.pdf  - figure: irregular looks
#   sims/futility-demo.pdf         - figure: futility trajectories
#   sims/recovery-eprocess.pdf     - figure: RECOVERY-like e-process
# =============================================================================

cat("=== Extended simulations for e-values paper ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Load evalinger package ---
# This script lives in evalinger/paper/; the package root is one level up.
pkg_dir <- if (file.exists(file.path("..", "DESCRIPTION"))) ".." else "."
cat("Loading evalinger from:", normalizePath(pkg_dir), "\n")
devtools::load_all(pkg_dir, quiet = TRUE)

out_dir <- "sims"
if (!dir.exists(out_dir)) dir.create(out_dir)

# =============================================================================
# 1. IRREGULAR LOOK SCHEDULE: GS VS E-VALUE
# =============================================================================
cat("--- 1. Irregular look schedule comparison ---\n")

# Design parameters (same as main study)
p_C     <- 0.30
p_T_alt <- 0.45
Nmax    <- 200
alpha   <- 0.025
nrep    <- 50000

lambda_star <- grow_lambda(p_T_alt, p_C)

# We compare three scenarios:
#   (a) Fixed 20 equally spaced looks (baseline)
#   (b) Irregular: 5 random looks drawn uniformly from [20, 200]
#   (c) Continuous: look at every patient pair (200 looks)
#
# For GS, we recalibrate the OBF boundary for each scenario.
# For e-values, the threshold stays at 1/alpha = 40 for all three.

set.seed(42)
xT_null <- matrix(rbinom(nrep * Nmax, 1, p_C),    nrow = nrep)
xC_null <- matrix(rbinom(nrep * Nmax, 1, p_C),    nrow = nrep)
xT_alt  <- matrix(rbinom(nrep * Nmax, 1, p_T_alt), nrow = nrep)
xC_alt  <- matrix(rbinom(nrep * Nmax, 1, p_C),    nrow = nrep)

# Cumulative log-e-values (computed once, used for all scenarios)
D_null <- xT_null - xC_null
D_alt  <- xT_alt  - xC_alt
log_inc_null <- log(pmax(1 + lambda_star * D_null, .Machine$double.xmin))
log_inc_alt  <- log(pmax(1 + lambda_star * D_alt,  .Machine$double.xmin))
cumlog_null <- t(apply(log_inc_null, 1, cumsum))
cumlog_alt  <- t(apply(log_inc_alt,  1, cumsum))

# Cumulative sums for z-statistics
cumT_null_full <- t(apply(xT_null, 1, cumsum))
cumC_null_full <- t(apply(xC_null, 1, cumsum))
cumT_alt_full  <- t(apply(xT_alt,  1, cumsum))
cumC_alt_full  <- t(apply(xC_alt,  1, cumsum))

threshold_log_e <- log(1 / alpha)

run_scenario <- function(look_times, label) {
  cat("  Scenario:", label, "(", length(look_times), "looks)\n")
  n_looks <- length(look_times)
  info_frac <- look_times / Nmax

  # Extract at look times
  cl_null <- cumlog_null[, look_times, drop = FALSE]
  cl_alt  <- cumlog_alt[, look_times,  drop = FALSE]
  cT_null <- cumT_null_full[, look_times, drop = FALSE]
  cC_null <- cumC_null_full[, look_times, drop = FALSE]
  cT_alt  <- cumT_alt_full[, look_times,  drop = FALSE]
  cC_alt  <- cumC_alt_full[, look_times,  drop = FALSE]
  nn <- matrix(look_times, nrow = nrep, ncol = n_looks, byrow = TRUE)

  # E-value: first crossing
  e_crossed_null <- cl_null >= threshold_log_e
  e_crossed_alt  <- cl_alt  >= threshold_log_e

  e_rej_null <- rowSums(e_crossed_null) > 0
  e_rej_alt  <- rowSums(e_crossed_alt)  > 0

  first_e_null <- apply(e_crossed_null, 1, function(r) {
    w <- which(r); if (length(w) > 0) look_times[w[1]] else Nmax
  })
  first_e_alt <- apply(e_crossed_alt, 1, function(r) {
    w <- which(r); if (length(w) > 0) look_times[w[1]] else Nmax
  })

  # GS: calibrate OBF boundary for this schedule
  z_null <- .z_matrix(cT_null, cC_null, nn)
  z_alt  <- .z_matrix(cT_alt,  cC_alt,  nn)

  scaled_z_null <- z_null * matrix(sqrt(info_frac), nrow = nrep,
                                    ncol = n_looks, byrow = TRUE)
  m_null <- apply(scaled_z_null, 1, max)
  gs_c <- as.numeric(quantile(m_null, probs = 1 - alpha, names = FALSE))
  obf_bounds <- gs_c / sqrt(info_frac)
  bounds_mat <- matrix(obf_bounds, nrow = nrep, ncol = n_looks, byrow = TRUE)

  gs_crossed_null <- z_null >= bounds_mat
  gs_crossed_alt  <- z_alt  >= bounds_mat

  gs_rej_null <- rowSums(gs_crossed_null) > 0
  gs_rej_alt  <- rowSums(gs_crossed_alt)  > 0

  first_gs_null <- apply(gs_crossed_null, 1, function(r) {
    w <- which(r); if (length(w) > 0) look_times[w[1]] else Nmax
  })
  first_gs_alt <- apply(gs_crossed_alt, 1, function(r) {
    w <- which(r); if (length(w) > 0) look_times[w[1]] else Nmax
  })

  data.frame(
    scenario = label,
    n_looks = n_looks,
    method = rep(c("E-value", "GS (recalibrated)"), each = 1),
    null_rej = c(mean(e_rej_null), mean(gs_rej_null)),
    alt_rej = c(mean(e_rej_alt), mean(gs_rej_alt)),
    avg_n_null = c(mean(first_e_null), mean(first_gs_null)),
    avg_n_alt = c(mean(first_e_alt), mean(first_gs_alt)),
    stringsAsFactors = FALSE
  )
}

# (a) Fixed 20 equally spaced
looks_fixed_20 <- round(seq(Nmax/20, Nmax, length.out = 20))
res_fixed_20 <- run_scenario(looks_fixed_20, "Fixed (20 looks)")

# (b) Fixed 5 equally spaced
looks_fixed_5 <- round(seq(Nmax/5, Nmax, length.out = 5))
res_fixed_5 <- run_scenario(looks_fixed_5, "Fixed (5 looks)")

# (c) Irregular: 5 random looks (repeated with 200 different random schedules,
#     then averaged to get expected performance under random scheduling)
cat("  Running 200 random schedules with 5 looks each...\n")
set.seed(123)
n_random_schedules <- 200
random_results <- vector("list", n_random_schedules)
for (s in seq_len(n_random_schedules)) {
  random_looks <- sort(c(sample(20:180, 4, replace = FALSE), Nmax))
  random_results[[s]] <- run_scenario(random_looks, "random")
}
# Average across random schedules
rand_df <- do.call(rbind, random_results)
res_irregular <- data.frame(
  scenario = "Irregular (5 random looks)",
  n_looks = 5,
  method = c("E-value", "GS (recalibrated)"),
  null_rej = tapply(rand_df$null_rej, rand_df$method, mean),
  alt_rej = tapply(rand_df$alt_rej, rand_df$method, mean),
  avg_n_null = tapply(rand_df$avg_n_null, rand_df$method, mean),
  avg_n_alt = tapply(rand_df$avg_n_alt, rand_df$method, mean),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# (d) Continuous: every patient pair
looks_continuous <- seq_len(Nmax)
res_continuous <- run_scenario(looks_continuous, "Continuous (200 looks)")

irregular_comparison <- rbind(res_fixed_5, res_fixed_20, res_irregular, res_continuous)
cat("  Results:\n")
print(irregular_comparison)

saveRDS(irregular_comparison, file.path(out_dir, "irregular_comparison.rds"))

# Figure: bar chart comparing power under different schedules
pdf(file.path(out_dir, "irregular-comparison.pdf"), width = 7, height = 4.5)
par(mar = c(5, 4.5, 2, 1))

# Extract power values
ev_power <- irregular_comparison$alt_rej[irregular_comparison$method == "E-value"]
gs_power <- irregular_comparison$alt_rej[irregular_comparison$method == "GS (recalibrated)"]
scenario_labels <- c("Fixed\n(5 looks)", "Fixed\n(20 looks)",
                     "Irregular\n(5 random)", "Continuous\n(200 looks)")

bp <- barplot(
  rbind(ev_power, gs_power),
  beside = TRUE,
  names.arg = scenario_labels,
  col = c("#D55E00", "#0072B2"),
  ylim = c(0, 1),
  ylab = "Power (rejection probability under H1)",
  main = expression("E-value vs group sequential power under different look schedules (" *
                    alpha == 0.025 * ")")
)
# Add Type I error as text
ev_t1e <- irregular_comparison$null_rej[irregular_comparison$method == "E-value"]
gs_t1e <- irregular_comparison$null_rej[irregular_comparison$method == "GS (recalibrated)"]
text(bp[1, ], ev_power + 0.03, sprintf("%.1f%%", 100 * ev_power), cex = 0.7)
text(bp[2, ], gs_power + 0.03, sprintf("%.1f%%", 100 * gs_power), cex = 0.7)

legend("topright",
       legend = c("E-value (fixed threshold)", "GS (recalibrated per schedule)"),
       fill = c("#D55E00", "#0072B2"),
       cex = 0.8, bty = "n")
dev.off()
cat("  Saved irregular-comparison.pdf\n")


# =============================================================================
# 2. PARAMETER GRID: VARY p_C, DELTA, N_LOOKS
# =============================================================================
cat("\n--- 2. Parameter grid ---\n")

grid_params <- expand.grid(
  p_C = c(0.10, 0.30, 0.50),
  delta = c(0.10, 0.15, 0.20),
  n_looks = c(5, 20),
  stringsAsFactors = FALSE
)
# Remove infeasible combinations where p_T >= 1
grid_params <- grid_params[grid_params$p_C + grid_params$delta < 1, ]

nrep_grid <- 20000
grid_results <- vector("list", nrow(grid_params))

for (i in seq_len(nrow(grid_params))) {
  pc <- grid_params$p_C[i]
  delta <- grid_params$delta[i]
  nl <- grid_params$n_looks[i]
  pt <- pc + delta
  cat(sprintf("  p_C=%.2f, delta=%.2f, n_looks=%d ... ", pc, delta, nl))

  cmp <- simulate_comparison(
    p_C = pc, p_T_alt = pt, Nmax = Nmax,
    n_looks = nl, alpha = alpha,
    methods = c("evalue", "gs_obf"),
    nrep = nrep_grid, seed = 42
  )

  lam <- grow_lambda(pt, pc)
  grate <- expected_growth_rate(lam, pt, pc)
  est_stop <- expected_stopping_time(lam, pt, pc, alpha)

  ev_row <- cmp$results[cmp$results$method == "evalue", ]
  gs_row <- cmp$results[cmp$results$method == "gs_obf", ]

  grid_results[[i]] <- data.frame(
    p_C = pc, delta = delta, p_T = pt,
    n_looks = nl,
    lambda_star = round(lam, 3),
    growth_rate = round(grate, 5),
    E_stop = round(est_stop),
    evalue_power = ev_row$alt_rej,
    evalue_t1e = ev_row$null_rej,
    gs_power = gs_row$alt_rej,
    gs_t1e = gs_row$null_rej,
    power_gap = gs_row$alt_rej - ev_row$alt_rej,
    stringsAsFactors = FALSE
  )
  cat(sprintf("E-val %.1f%% vs GS %.1f%% (gap %.1f pp)\n",
              100 * ev_row$alt_rej, 100 * gs_row$alt_rej,
              100 * (gs_row$alt_rej - ev_row$alt_rej)))
}

parameter_grid <- do.call(rbind, grid_results)
cat("  Parameter grid results:\n")
print(parameter_grid)
saveRDS(parameter_grid, file.path(out_dir, "parameter_grid.rds"))


# =============================================================================
# 3. FUTILITY MONITORING DEMONSTRATION
# =============================================================================
cat("\n--- 3. Futility monitoring demonstration ---\n")

# Scenario: treatment is INEFFECTIVE (small effect, below MCID)
# p_T = 0.33, p_C = 0.30, delta = 0.03 (true), MCID = 0.10
p_C_fut <- 0.30
p_T_fut <- 0.33   # true effect well below MCID
delta_min <- 0.10  # minimum clinically important difference
Nmax_fut <- 300

# Run nrep simulations to get futility detection rates
nrep_fut <- 10000
set.seed(2024)

fut_cs_detected <- integer(nrep_fut)
fut_ep_detected <- integer(nrep_fut)
lambda_fut <- grow_lambda(p_C_fut + delta_min, p_C_fut)  # tuned for MCID

for (rep in seq_len(nrep_fut)) {
  x_T <- rbinom(Nmax_fut, 1, p_T_fut)
  x_C <- rbinom(Nmax_fut, 1, p_C_fut)

  # CS-based futility
  cs <- confseq_binary(x_T, x_C, alpha = 0.05, t_min = 20)
  fc <- futility_cs(cs, delta_min = delta_min)
  fut_cs_detected[rep] <- ifelse(is.na(fc$first_futile), Nmax_fut, fc$first_futile)

  # Reciprocal e-process futility
  fe <- futility_eprocess(x_T, x_C, delta_min = delta_min,
                          alpha_f = 0.10, lambda_f = 0.3)
  fut_ep_detected[rep] <- ifelse(is.na(fe$first_futile), Nmax_fut, fe$first_futile)
}

futility_demo <- list(
  design = list(p_T = p_T_fut, p_C = p_C_fut, delta_min = delta_min,
                Nmax = Nmax_fut, nrep = nrep_fut),
  cs_futility_rate = mean(fut_cs_detected < Nmax_fut),
  cs_median_time = median(fut_cs_detected[fut_cs_detected < Nmax_fut]),
  cs_mean_time = mean(fut_cs_detected),
  ep_futility_rate = mean(fut_ep_detected < Nmax_fut),
  ep_median_time = median(fut_ep_detected[fut_ep_detected < Nmax_fut]),
  ep_mean_time = mean(fut_ep_detected)
)

cat(sprintf("  CS-based futility: detected %.1f%%, median time %.0f (among detected)\n",
            100 * futility_demo$cs_futility_rate,
            futility_demo$cs_median_time))
cat(sprintf("  E-process futility: detected %.1f%%, median time %.0f (among detected)\n",
            100 * futility_demo$ep_futility_rate,
            futility_demo$ep_median_time))

saveRDS(futility_demo, file.path(out_dir, "futility_demo.rds"))

# Figure: single-trial futility illustration
set.seed(77)
x_T_ex <- rbinom(Nmax_fut, 1, p_T_fut)
x_C_ex <- rbinom(Nmax_fut, 1, p_C_fut)
cs_ex <- confseq_binary(x_T_ex, x_C_ex, alpha = 0.05, t_min = 20)
fc_ex <- futility_cs(cs_ex, delta_min = delta_min)
fe_ex <- futility_eprocess(x_T_ex, x_C_ex, delta_min = delta_min,
                           alpha_f = 0.10, lambda_f = 0.3)

pdf(file.path(out_dir, "futility-demo.pdf"), width = 9, height = 4)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))

# Panel 1: Confidence sequence with MCID
idx <- 20:Nmax_fut
plot(idx, cs_ex$delta_hat[idx], type = "l", lwd = 2,
     ylim = range(c(cs_ex$lower[idx], cs_ex$upper[idx], delta_min)),
     xlab = "Patient pairs", ylab = expression(hat(delta) == hat(p)[T] - hat(p)[C]),
     main = "Confidence sequence futility")
polygon(c(idx, rev(idx)),
        c(cs_ex$lower[idx], rev(cs_ex$upper[idx])),
        col = adjustcolor("steelblue", 0.2), border = NA)
lines(idx, cs_ex$lower[idx], col = "steelblue", lty = 2)
lines(idx, cs_ex$upper[idx], col = "steelblue", lty = 2)
lines(idx, cs_ex$delta_hat[idx], lwd = 2)
abline(h = delta_min, col = "red", lty = 3, lwd = 1.5)
abline(h = 0, col = "gray80")
if (!is.na(fc_ex$first_futile)) {
  abline(v = fc_ex$first_futile, col = "red", lty = 2)
  text(fc_ex$first_futile + 5, delta_min + 0.02,
       paste0("Futility at n=", fc_ex$first_futile),
       adj = 0, cex = 0.7, col = "red")
}
legend("topright",
       legend = c("Point estimate", "95% CS",
                  expression(delta[min] == 0.10)),
       col = c("black", "steelblue", "red"),
       lwd = c(2, 1, 1.5), lty = c(1, 2, 3), cex = 0.7, bty = "n")

# Panel 2: Reciprocal e-process
plot(seq_len(Nmax_fut), fe_ex$log_evalue, type = "l", lwd = 2,
     col = "#D55E00",
     xlab = "Patient pairs", ylab = expression(log ~ E[n]^"'"),
     main = "Reciprocal e-process futility")
abline(h = log(1 / 0.10), lty = 3, col = "black", lwd = 1.5)
abline(h = 0, col = "gray80")
if (!is.na(fe_ex$first_futile)) {
  abline(v = fe_ex$first_futile, col = "red", lty = 2)
  text(fe_ex$first_futile + 5, log(1/0.10) * 0.5,
       paste0("Futility at n=", fe_ex$first_futile),
       adj = 0, cex = 0.7, col = "red")
}
legend("topleft",
       legend = c("Reciprocal e-process",
                  expression(paste("Threshold ", log(1/alpha[f])))),
       col = c("#D55E00", "black"),
       lwd = c(2, 1.5), lty = c(1, 3), cex = 0.7, bty = "n")
dev.off()
cat("  Saved futility-demo.pdf\n")


# =============================================================================
# 4. RECOVERY-LIKE LARGE-TRIAL SIMULATION
# =============================================================================
cat("\n--- 4. RECOVERY-like simulation ---\n")
# Published RECOVERY trial data (dexamethasone vs usual care, 28-day mortality):
#   Dexamethasone: 482/2104 deaths = 22.9%
#   Usual care:   1110/4321 deaths = 25.7%
# We simulate a simplified version with 1:1 randomization
# using the observed mortality rates.

p_C_rec <- 0.257   # usual care mortality
p_T_rec <- 0.229   # dexamethasone mortality
N_rec   <- 2000     # per arm (simplified 1:1 from the ~6400 total)
delta_rec <- p_C_rec - p_T_rec  # ~2.8pp mortality reduction

# Design: test H_0: p_T >= p_C vs H_1: p_T < p_C (treatment reduces mortality)
# Equivalently, test p_C > p_T using x_C as "treatment" and x_T as "control"
# in the eprocess_binary framework (which tests first arg > second arg)
lambda_rec <- grow_lambda(p_C_rec, p_T_rec)
g_rec <- expected_growth_rate(lambda_rec, p_C_rec, p_T_rec)
est_stop_rec <- expected_stopping_time(lambda_rec, p_C_rec, p_T_rec, alpha = alpha)

cat(sprintf("  RECOVERY design: p_UC=%.3f, p_Dex=%.3f, delta=%.3f\n",
            p_C_rec, p_T_rec, delta_rec))
cat(sprintf("  GROW lambda: %.4f, growth rate: %.6f nats/pair\n",
            lambda_rec, g_rec))
cat(sprintf("  Expected stopping time: %.0f pairs\n", est_stop_rec))

# Simulate a single large RECOVERY-like trial
set.seed(2020)  # year of RECOVERY
x_UC  <- rbinom(N_rec, 1, p_C_rec)  # usual care (higher mortality)
x_Dex <- rbinom(N_rec, 1, p_T_rec)  # dexamethasone (lower mortality)

# E-process: test UC > Dex (mortality reduction)
ep_rec <- eprocess_binary(x_UC, x_Dex, lambda = lambda_rec, alpha = alpha)

# Also run with monthly looks (every ~167 patients, 12 looks)
look_rec <- round(seq(N_rec/12, N_rec, length.out = 12))
hm_rec <- hybrid_monitor(x_UC, x_Dex, look_times = look_rec,
                         lambda = lambda_rec, alpha = alpha,
                         p_T_design = p_C_rec, p_C_design = p_T_rec)

# Confidence sequence
cs_rec <- confseq_binary(x_UC, x_Dex, alpha = 0.05, t_min = 50)

# Monte Carlo: rejection rate and stopping time
nrep_rec <- 10000
rec_rej <- logical(nrep_rec)
rec_stop <- integer(nrep_rec)
set.seed(2020)
for (rep in seq_len(nrep_rec)) {
  xuc <- rbinom(N_rec, 1, p_C_rec)
  xdx <- rbinom(N_rec, 1, p_T_rec)
  D <- xuc - xdx
  log_e <- cumsum(log(pmax(1 + lambda_rec * D, .Machine$double.xmin)))
  crossed <- log_e >= threshold_log_e
  if (any(crossed)) {
    rec_rej[rep] <- TRUE
    rec_stop[rep] <- which(crossed)[1]
  } else {
    rec_stop[rep] <- N_rec
  }
}

recovery_sim <- list(
  design = list(p_UC = p_C_rec, p_Dex = p_T_rec, N_per_arm = N_rec,
                lambda = lambda_rec, growth_rate = g_rec,
                est_stop = est_stop_rec),
  single_trial = list(
    rejected = ep_rec$rejected,
    rejection_time = ep_rec$rejection_time,
    final_log_e = ep_rec$log_evalue[N_rec],
    final_e = ep_rec$evalue[N_rec],
    cs_final = list(
      estimate = cs_rec$delta_hat[N_rec],
      lower = cs_rec$lower[N_rec],
      upper = cs_rec$upper[N_rec]
    )
  ),
  hybrid = hm_rec,
  monte_carlo = list(
    nrep = nrep_rec,
    power = mean(rec_rej),
    avg_stop_alt = mean(rec_stop),
    median_stop_alt = median(rec_stop[rec_rej])
  )
)

cat(sprintf("  Single trial: E = %.1f, rejected = %s at n = %s\n",
            ep_rec$evalue[N_rec],
            ep_rec$rejected,
            ifelse(is.na(ep_rec$rejection_time), "NA", ep_rec$rejection_time)))
cat(sprintf("  CS at n=%d: [%.4f, %.4f], est = %.4f\n",
            N_rec, cs_rec$lower[N_rec], cs_rec$upper[N_rec], cs_rec$delta_hat[N_rec]))
cat(sprintf("  Monte Carlo power (n=%d per arm): %.1f%%\n",
            N_rec, 100 * mean(rec_rej)))
cat(sprintf("  Median stopping time (among rejections): %d\n",
            median(rec_stop[rec_rej])))

saveRDS(recovery_sim, file.path(out_dir, "recovery_sim.rds"))

# Figure: e-process trajectory for RECOVERY-like trial
pdf(file.path(out_dir, "recovery-eprocess.pdf"), width = 8, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))

# Panel 1: E-process trajectory
nn <- seq_len(N_rec)
plot(nn, ep_rec$log_evalue, type = "l", lwd = 1.5,
     col = "#D55E00",
     xlab = "Patient pairs",
     ylab = expression(log ~ E[n]),
     main = "RECOVERY-like trial: e-process")
abline(h = threshold_log_e, lty = 3, col = "black", lwd = 1.5)
abline(h = 0, col = "gray80")
if (ep_rec$rejected) {
  abline(v = ep_rec$rejection_time, col = "red", lty = 2)
  text(ep_rec$rejection_time, threshold_log_e * 0.6,
       paste0("Reject at n=", ep_rec$rejection_time),
       adj = c(-0.1, 0), cex = 0.7, col = "red")
}
legend("topleft",
       legend = c(
         sprintf("E-process (%.1f%% vs %.1f%%)", 100*p_C_rec, 100*p_T_rec),
         expression(paste("Threshold  ", log(1/alpha)))
       ),
       col = c("#D55E00", "black"),
       lwd = c(1.5, 1.5), lty = c(1, 3), cex = 0.7, bty = "n")

# Panel 2: Confidence sequence
idx2 <- 50:N_rec
plot(idx2, cs_rec$delta_hat[idx2], type = "l", lwd = 2,
     ylim = range(c(cs_rec$lower[idx2], cs_rec$upper[idx2])),
     xlab = "Patient pairs",
     ylab = expression(hat(delta) == hat(p)[UC] - hat(p)[Dex]),
     main = "RECOVERY-like trial: confidence sequence")
polygon(c(idx2, rev(idx2)),
        c(cs_rec$lower[idx2], rev(cs_rec$upper[idx2])),
        col = adjustcolor("steelblue", 0.2), border = NA)
lines(idx2, cs_rec$lower[idx2], col = "steelblue", lty = 2)
lines(idx2, cs_rec$upper[idx2], col = "steelblue", lty = 2)
lines(idx2, cs_rec$delta_hat[idx2], lwd = 2)
abline(h = 0, lty = 3, col = "gray50")
abline(h = delta_rec, lty = 3, col = "darkgreen")
legend("topright",
       legend = c("Point estimate", "95% CS",
                  sprintf("True delta = %.3f", delta_rec)),
       col = c("black", "steelblue", "darkgreen"),
       lwd = c(2, 1, 1), lty = c(1, 2, 3), cex = 0.7, bty = "n")
dev.off()
cat("  Saved recovery-eprocess.pdf\n")


# =============================================================================
cat("\n=== All extended simulations complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("Outputs saved to:", normalizePath(out_dir), "\n")
