#' Simulate and Compare Monitoring Methods
#'
#' Run a Monte Carlo simulation comparing e-value monitoring with group
#' sequential, naive repeated testing, and Bayesian monitoring approaches
#' for two-arm binary trials.
#'
#' @param p_C Control arm response probability.
#' @param p_T_alt Treatment arm response probability under the alternative.
#' @param Nmax Maximum sample size per arm.
#' @param n_looks Number of equally-spaced interim looks.
#' @param alpha Significance level (default 0.025).
#' @param methods Character vector of methods to compare. Options:
#'   \code{"evalue"}, \code{"gs_obf"} (O'Brien-Fleming-like group sequential),
#'   \code{"naive_p"} (repeated unadjusted p-values),
#'   \code{"cal_bayes"} (calibrated Bayesian posterior probability).
#'   Default includes all four.
#' @param nrep Number of Monte Carlo replications (default 5000).
#' @param seed Random seed (default 42).
#' @param lambda Optional fixed betting fraction for the e-value method. If
#'   NULL, uses the GROW-optimal fraction.
#'
#' @return An object of class \code{"ecomparison"} with components:
#' \describe{
#'   \item{results}{A data frame with columns \code{method}, \code{null_rej},
#'     \code{alt_rej}, \code{avg_n_null}, \code{avg_n_alt}.}
#'   \item{design}{A list with the simulation parameters.}
#'   \item{raw}{A list of per-method raw results for further analysis.}
#' }
#'
#' @details
#' Under the null, both arms draw from \code{Bernoulli(p_C)}. Under the
#' alternative, the treatment arm draws from \code{Bernoulli(p_T_alt)}.
#'
#' The group sequential method uses an O'Brien-Fleming-like alpha spending
#' function with equally-spaced looks. The calibrated Bayesian method uses
#' a \code{Beta(0.5, 0.5)} (Jeffreys) prior on each arm, computes the
#' posterior probability \eqn{P(p_T > p_C \mid \text{data})} analytically,
#' and calibrates the decision threshold to achieve the nominal Type I error.
#'
#' @examples
#' cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
#'                            n_looks = 4, nrep = 1000, seed = 1)
#' cmp
#'
#' @export
simulate_comparison <- function(p_C, p_T_alt, Nmax, n_looks = 4,
                                alpha = 0.025,
                                methods = c("evalue", "gs_obf", "naive_p",
                                            "cal_bayes"),
                                nrep = 5000, seed = 42, lambda = NULL) {
  stopifnot(p_T_alt > p_C, p_T_alt < 1, p_C > 0, p_C < 1)
  stopifnot(Nmax >= n_looks * 2)
  stopifnot(n_looks >= 1)
  stopifnot(alpha > 0, alpha < 1)

  if (is.null(lambda)) {
    lambda <- grow_lambda(p_T_alt, p_C)
  }

  look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))
  threshold_log_e <- log(1 / alpha)
  obf_bounds <- .obf_boundaries(look_times, Nmax, alpha)

  # --- Generate all trial data upfront as matrices (nrep x Nmax) ---
  set.seed(seed)
  xT_null_mat <- matrix(stats::rbinom(nrep * Nmax, 1, p_C),    nrow = nrep)
  xC_null_mat <- matrix(stats::rbinom(nrep * Nmax, 1, p_C),    nrow = nrep)
  xT_alt_mat  <- matrix(stats::rbinom(nrep * Nmax, 1, p_T_alt), nrow = nrep)
  xC_alt_mat  <- matrix(stats::rbinom(nrep * Nmax, 1, p_C),    nrow = nrep)

  # --- Pre-compute cumulative statistics used by multiple methods ---
  # Cumulative sums at look_times (nrep x n_looks matrices)
  cumT_null <- t(apply(xT_null_mat, 1, cumsum))[, look_times, drop = FALSE]
  cumC_null <- t(apply(xC_null_mat, 1, cumsum))[, look_times, drop = FALSE]
  cumT_alt  <- t(apply(xT_alt_mat,  1, cumsum))[, look_times, drop = FALSE]
  cumC_alt  <- t(apply(xC_alt_mat,  1, cumsum))[, look_times, drop = FALSE]

  # Cumulative log-e-values for evalue method (nrep x Nmax matrices)
  if ("evalue" %in% methods) {
    D_null <- xT_null_mat - xC_null_mat
    D_alt  <- xT_alt_mat  - xC_alt_mat
    log_inc_null <- log(pmax(1 + lambda * D_null, .Machine$double.xmin))
    log_inc_alt  <- log(pmax(1 + lambda * D_alt,  .Machine$double.xmin))
    cum_log_e_null <- t(apply(log_inc_null, 1, cumsum))[, look_times, drop = FALSE]
    cum_log_e_alt  <- t(apply(log_inc_alt,  1, cumsum))[, look_times, drop = FALSE]
  }

  # --- Calibrate Bayesian threshold (if needed) ---
  if ("cal_bayes" %in% methods) {
    bayes_thresh <- .calibrate_bayes_threshold_vec(
      cumT_null, cumC_null, look_times, alpha,
      nrep_cal = min(nrep, 5000)
    )
  }

  # --- Run each method in vectorized fashion ---
  raw <- list()

  for (m in methods) {
    rej_null <- logical(nrep)
    rej_alt  <- logical(nrep)
    stop_null <- rep(Nmax, nrep)
    stop_alt  <- rep(Nmax, nrep)

    for (scenario in c("null", "alt")) {
      if (scenario == "null") {
        cumT <- cumT_null; cumC <- cumC_null
      } else {
        cumT <- cumT_alt; cumC <- cumC_alt
      }

      # For each rep, find first look that crosses the threshold
      if (m == "evalue") {
        cum_le <- if (scenario == "null") cum_log_e_null else cum_log_e_alt
        crossed <- cum_le >= threshold_log_e  # nrep x n_looks logical matrix
        .store_results(crossed, look_times, Nmax, scenario, nrep,
                       rej_null, rej_alt, stop_null, stop_alt) -> tmp
        rej_null <- tmp$rej_null; rej_alt <- tmp$rej_alt
        stop_null <- tmp$stop_null; stop_alt <- tmp$stop_alt

      } else if (m == "gs_obf") {
        # Z-statistics at each look: delta_hat / se
        nn <- matrix(look_times, nrow = nrep, ncol = length(look_times), byrow = TRUE)
        phat_T <- cumT / nn
        phat_C <- cumC / nn
        delta <- phat_T - phat_C
        se <- sqrt(phat_T * (1 - phat_T) / nn + phat_C * (1 - phat_C) / nn)
        se[se < 1e-12] <- 1e-6
        z <- delta / se
        bounds_mat <- matrix(obf_bounds, nrow = nrep, ncol = length(look_times), byrow = TRUE)
        crossed <- z >= bounds_mat
        .store_results(crossed, look_times, Nmax, scenario, nrep,
                       rej_null, rej_alt, stop_null, stop_alt) -> tmp
        rej_null <- tmp$rej_null; rej_alt <- tmp$rej_alt
        stop_null <- tmp$stop_null; stop_alt <- tmp$stop_alt

      } else if (m == "naive_p") {
        nn <- matrix(look_times, nrow = nrep, ncol = length(look_times), byrow = TRUE)
        phat_T <- cumT / nn
        phat_C <- cumC / nn
        delta <- phat_T - phat_C
        se <- sqrt(phat_T * (1 - phat_T) / nn + phat_C * (1 - phat_C) / nn)
        se[se < 1e-12] <- 1e-6
        z <- delta / se
        pvals <- stats::pnorm(z, lower.tail = FALSE)
        crossed <- pvals < alpha
        .store_results(crossed, look_times, Nmax, scenario, nrep,
                       rej_null, rej_alt, stop_null, stop_alt) -> tmp
        rej_null <- tmp$rej_null; rej_alt <- tmp$rej_alt
        stop_null <- tmp$stop_null; stop_alt <- tmp$stop_alt

      } else if (m == "cal_bayes") {
        nn <- matrix(look_times, nrow = nrep, ncol = length(look_times), byrow = TRUE)
        # Analytic P(p_T > p_C | data) for independent Beta posteriors
        pp <- .posterior_prob_matrix(cumT, cumC, nn)
        crossed <- pp >= bayes_thresh
        .store_results(crossed, look_times, Nmax, scenario, nrep,
                       rej_null, rej_alt, stop_null, stop_alt) -> tmp
        rej_null <- tmp$rej_null; rej_alt <- tmp$rej_alt
        stop_null <- tmp$stop_null; stop_alt <- tmp$stop_alt
      }
    }

    raw[[m]] <- list(
      null_rej = rej_null, alt_rej = rej_alt,
      null_stop = stop_null, alt_stop = stop_alt
    )
  }

  # --- Summarize ---
  rows <- lapply(methods, function(m) {
    data.frame(
      method = m,
      null_rej = mean(raw[[m]]$null_rej),
      alt_rej  = mean(raw[[m]]$alt_rej),
      avg_n_null = mean(raw[[m]]$null_stop),
      avg_n_alt  = mean(raw[[m]]$alt_stop),
      stringsAsFactors = FALSE
    )
  })
  results <- do.call(rbind, rows)

  structure(
    list(
      results = results,
      design = list(p_C = p_C, p_T_alt = p_T_alt, Nmax = Nmax,
                    n_looks = n_looks, alpha = alpha, nrep = nrep,
                    lambda = lambda),
      raw = raw
    ),
    class = "ecomparison"
  )
}

# --- Internal helpers ---

#' Extract first crossing time from a logical matrix (nrep x n_looks)
#' @keywords internal
.store_results <- function(crossed, look_times, Nmax, scenario, nrep,
                           rej_null, rej_alt, stop_null, stop_alt) {
  # For each row, find first TRUE column (first crossing)
  any_crossed <- rowSums(crossed) > 0
  # first_cross_col: column index of first TRUE per row (NA if none)
  first_col <- max.col(crossed, ties.method = "first")
  first_col[!any_crossed] <- NA_integer_
  stop_time <- ifelse(any_crossed, look_times[first_col], Nmax)

  if (scenario == "null") {
    rej_null <- any_crossed
    stop_null <- stop_time
  } else {
    rej_alt <- any_crossed
    stop_alt <- stop_time
  }
  list(rej_null = rej_null, rej_alt = rej_alt,
       stop_null = stop_null, stop_alt = stop_alt)
}

#' Analytic P(p_T > p_C | data) for independent Beta posteriors (vectorized)
#'
#' Uses the normal approximation to the Beta posterior difference, which is
#' fast and accurate for n >= 10. For small n, the approximation is slightly
#' less accurate but sufficient for simulation-based calibration.
#' @keywords internal
.posterior_prob_matrix <- function(cumT, cumC, nn) {
  # Jeffreys prior: Beta(0.5, 0.5)
  a_T <- 0.5 + cumT
  b_T <- 0.5 + nn - cumT
  a_C <- 0.5 + cumC
  b_C <- 0.5 + nn - cumC

  # Normal approximation to Beta posteriors:
  # E[p_T] = a_T/(a_T+b_T), Var[p_T] = a_T*b_T / ((a_T+b_T)^2 * (a_T+b_T+1))
  n_T <- a_T + b_T
  n_C <- a_C + b_C
  mu_T <- a_T / n_T
  mu_C <- a_C / n_C
  var_T <- a_T * b_T / (n_T^2 * (n_T + 1))
  var_C <- a_C * b_C / (n_C^2 * (n_C + 1))

  # P(p_T > p_C) = P(p_T - p_C > 0) ~ Phi(mu_diff / sd_diff)
  mu_diff <- mu_T - mu_C
  sd_diff <- sqrt(var_T + var_C)
  sd_diff[sd_diff < 1e-12] <- 1e-12
  stats::pnorm(mu_diff / sd_diff)
}

#' Calibrate Bayesian threshold using vectorized approach
#' @keywords internal
.calibrate_bayes_threshold_vec <- function(cumT_null, cumC_null, look_times,
                                           alpha, nrep_cal = 5000) {
  # Use at most nrep_cal rows from the null data
  nr <- min(nrow(cumT_null), nrep_cal)
  cumT <- cumT_null[seq_len(nr), , drop = FALSE]
  cumC <- cumC_null[seq_len(nr), , drop = FALSE]
  nn <- matrix(look_times, nrow = nr, ncol = length(look_times), byrow = TRUE)

  pp <- .posterior_prob_matrix(cumT, cumC, nn)
  max_pp <- apply(pp, 1, max)
  stats::quantile(max_pp, 1 - alpha, names = FALSE)
}

#' @keywords internal
.obf_boundaries <- function(look_times, Nmax, alpha) {
  info_frac <- look_times / Nmax
  z_alpha <- stats::qnorm(1 - alpha)
  z_alpha / sqrt(info_frac)
}

#' @keywords internal
.posterior_prob_greater <- function(a_T, b_T, a_C, b_C, nsim = 5000) {
  # Keep for backward compatibility / standalone use
  n_T <- a_T + b_T
  n_C <- a_C + b_C
  mu_T <- a_T / n_T
  mu_C <- a_C / n_C
  var_T <- a_T * b_T / (n_T^2 * (n_T + 1))
  var_C <- a_C * b_C / (n_C^2 * (n_C + 1))
  mu_diff <- mu_T - mu_C
  sd_diff <- sqrt(var_T + var_C)
  if (sd_diff < 1e-12) sd_diff <- 1e-12
  stats::pnorm(mu_diff / sd_diff)
}

# --- Legacy wrappers for backward compatibility ---
#' @keywords internal
.run_method <- function(method, x_T, x_C, look_times, lambda,
                        threshold_log_e, obf_bounds, bayes_thresh,
                        alpha, Nmax) {
  switch(method,
    evalue = .run_evalue(x_T, x_C, look_times, lambda, threshold_log_e, Nmax),
    gs_obf = .run_gs(x_T, x_C, look_times, obf_bounds, Nmax),
    naive_p = .run_naive_p(x_T, x_C, look_times, alpha, Nmax),
    cal_bayes = .run_bayes(x_T, x_C, look_times, bayes_thresh, Nmax),
    stop("Unknown method: ", method)
  )
}

#' @keywords internal
.run_evalue <- function(x_T, x_C, look_times, lambda, threshold, Nmax) {
  D <- as.numeric(x_T) - as.numeric(x_C)
  log_e <- cumsum(log(pmax(1 + lambda * D, .Machine$double.xmin)))
  for (t in look_times) {
    if (log_e[t] >= threshold) return(list(rejected = TRUE, stop_n = t))
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.run_gs <- function(x_T, x_C, look_times, obf_bounds, Nmax) {
  for (k in seq_along(look_times)) {
    t <- look_times[k]
    p_hat_T <- sum(x_T[1:t]) / t
    p_hat_C <- sum(x_C[1:t]) / t
    se <- sqrt(p_hat_T * (1 - p_hat_T) / t + p_hat_C * (1 - p_hat_C) / t)
    if (se < 1e-12) se <- 1e-6
    z <- (p_hat_T - p_hat_C) / se
    if (z >= obf_bounds[k]) return(list(rejected = TRUE, stop_n = t))
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.run_naive_p <- function(x_T, x_C, look_times, alpha, Nmax) {
  for (t in look_times) {
    p_hat_T <- sum(x_T[1:t]) / t
    p_hat_C <- sum(x_C[1:t]) / t
    se <- sqrt(p_hat_T * (1 - p_hat_T) / t + p_hat_C * (1 - p_hat_C) / t)
    if (se < 1e-12) se <- 1e-6
    z <- (p_hat_T - p_hat_C) / se
    pval <- stats::pnorm(z, lower.tail = FALSE)
    if (pval < alpha) return(list(rejected = TRUE, stop_n = t))
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.run_bayes <- function(x_T, x_C, look_times, threshold, Nmax) {
  for (t in look_times) {
    n_T <- sum(x_T[1:t])
    n_C <- sum(x_C[1:t])
    a_T <- 0.5 + n_T; b_T <- 0.5 + t - n_T
    a_C <- 0.5 + n_C; b_C <- 0.5 + t - n_C
    pp <- .posterior_prob_greater(a_T, b_T, a_C, b_C)
    if (pp >= threshold) return(list(rejected = TRUE, stop_n = t))
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.calibrate_bayes_threshold <- function(p_C, look_times, alpha,
                                       nrep_cal = 5000, seed = 999) {
  set.seed(seed)
  Nmax <- max(look_times)
  max_pp <- numeric(nrep_cal)
  for (rep in seq_len(nrep_cal)) {
    x_T <- stats::rbinom(Nmax, 1, p_C)
    x_C <- stats::rbinom(Nmax, 1, p_C)
    pp_vec <- numeric(length(look_times))
    for (k in seq_along(look_times)) {
      t <- look_times[k]
      n_T <- sum(x_T[1:t]); n_C <- sum(x_C[1:t])
      a_T <- 0.5 + n_T; b_T <- 0.5 + t - n_T
      a_C <- 0.5 + n_C; b_C <- 0.5 + t - n_C
      pp_vec[k] <- .posterior_prob_greater(a_T, b_T, a_C, b_C)
    }
    max_pp[rep] <- max(pp_vec)
  }
  stats::quantile(max_pp, 1 - alpha, names = FALSE)
}
