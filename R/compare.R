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
#'   \code{"naive_bayes"} (naive Bayesian posterior threshold),
#'   \code{"cal_bayes"} (calibrated Bayesian posterior probability).
#'   Default includes all five.
#' @param nrep Number of Monte Carlo replications (default 5000).
#' @param seed Random seed (default 42).
#' @param lambda Optional fixed betting fraction for the e-value method. If
#'   NULL, uses the GROW-optimal fraction.
#' @param posterior_method Method for computing the Bayesian posterior
#'   probability \eqn{P(p_T > p_C \mid \text{data})} when
#'   \code{"cal_bayes"} is included in \code{methods}. One of:
#'   \describe{
#'     \item{\code{"normal"}}{(Default.) Normal approximation to the
#'       difference of independent Beta posteriors. Each Beta posterior is
#'       approximated by its mean and variance, and the probability is
#'       obtained via \code{\link[stats]{pnorm}}.
#'       Highly accurate for per-arm sample sizes \eqn{\ge 10}; extremely fast
#'       because it operates on entire matrices without any loops or numerical
#'       integration.}
#'     \item{\code{"exact"}}{Exact computation using the regularized
#'       incomplete beta function. For independent posteriors
#'       \eqn{p_T \sim \text{Beta}(a_T, b_T)} and
#'       \eqn{p_C \sim \text{Beta}(a_C, b_C)}, the probability
#'       \deqn{P(p_T > p_C) = \int_0^1 I_x(a_C, b_C)\, f_{a_T,b_T}(x)\,dx}
#'       where \eqn{I_x(a,b)} is the regularized incomplete beta function
#'       (\code{\link[stats]{pbeta}}) and \eqn{f_{a_T,b_T}} is the Beta
#'       density.
#'       This is computed via \code{\link[stats]{integrate}} and is exact to
#'       machine precision, but substantially slower than the normal
#'       approximation because each element requires scalar numerical
#'       integration.}
#'   }
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
#' The group sequential method uses an O'Brien--Fleming-like boundary of the
#' form \eqn{z \ge c/\sqrt{t}} (where \eqn{t} is information fraction) and
#' calibrates the constant \eqn{c} by simulation to achieve the nominal Type I
#' error across the \code{n_looks} interim looks. The calibrated Bayesian method
#' uses a \code{Beta(0.5, 0.5)} (Jeffreys) prior on each arm, computes the
#' posterior probability \eqn{P(p_T > p_C \mid \text{data})}, and calibrates the
#' posterior threshold to achieve the nominal Type I error across looks.
#'
#' \strong{Posterior probability computation.}
#' The posterior probability \eqn{P(p_T > p_C \mid \text{data})} for
#' independent Beta posteriors can be computed in two ways (controlled by
#' \code{posterior_method}):
#'
#' \enumerate{
#'   \item \strong{Normal approximation} (\code{"normal"}, default): Each Beta
#'     posterior is approximated by a Gaussian with matched mean and variance.
#'     The probability reduces to a single \code{pnorm} call per cell. This is
#'     extremely fast and suitable for the simulation setting where per-arm
#'     samples are typically \eqn{\ge 10}.
#'
#'   \item \strong{Exact via regularized incomplete beta function}
#'     (\code{"exact"}): The exact probability is obtained by numerically
#'     integrating \eqn{I_x(a_C, b_C) \cdot f_{a_T, b_T}(x)} over
#'     \eqn{[0, 1]}, where \eqn{I_x(a, b) = } \code{pbeta(x, a, b)} is the
#'     regularized incomplete beta function. This is exact to machine
#'     precision but requires scalar \code{integrate} calls, making it
#'     substantially slower for large-scale simulations.
#' }
#'
#' For most simulation purposes, \code{"normal"} is recommended. Use
#' \code{"exact"} when you need machine-precision posterior probabilities,
#' e.g., for validating results or when sample sizes per look are very small.
#'
#' @examples
#' cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
#'                            n_looks = 4, nrep = 1000, seed = 1)
#' cmp
#'
#' # Use exact posterior probability via regularized incomplete beta function
#' cmp_exact <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
#'                                  n_looks = 4, nrep = 500, seed = 1,
#'                                  posterior_method = "exact")
#'
#' @export
simulate_comparison <- function(p_C, p_T_alt, Nmax, n_looks = 4,
                                alpha = 0.025,
                                methods = c("evalue", "gs_obf", "naive_p",
                                            "naive_bayes", "cal_bayes"),
                                nrep = 5000, seed = 42, lambda = NULL,
                                posterior_method = c("normal", "exact")) {
  posterior_method <- match.arg(posterior_method)
  stopifnot(p_T_alt > p_C, p_T_alt < 1, p_C > 0, p_C < 1)
  stopifnot(Nmax >= n_looks * 2)
  stopifnot(n_looks >= 1)
  stopifnot(alpha > 0, alpha < 1)

  if (is.null(lambda)) {
    lambda <- grow_lambda(p_T_alt, p_C)
  }

  look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))
  threshold_log_e <- log(1 / alpha)
  info_frac <- look_times / Nmax
  obf_bounds <- NULL
  gs_c <- NULL

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

  # Shared nn matrix for interim looks
  nn_looks <- matrix(look_times, nrow = nrep, ncol = length(look_times), byrow = TRUE)

  # Cumulative log-e-values for evalue method (nrep x Nmax matrices)
  if ("evalue" %in% methods) {
    D_null <- xT_null_mat - xC_null_mat
    D_alt  <- xT_alt_mat  - xC_alt_mat
    log_inc_null <- log(pmax(1 + lambda * D_null, .Machine$double.xmin))
    log_inc_alt  <- log(pmax(1 + lambda * D_alt,  .Machine$double.xmin))
    cum_log_e_null <- t(apply(log_inc_null, 1, cumsum))[, look_times, drop = FALSE]
    cum_log_e_alt  <- t(apply(log_inc_alt,  1, cumsum))[, look_times, drop = FALSE]
  }

  # Pre-compute z-statistics for methods using Wald z (GS and naive p)
  if (any(c("gs_obf", "naive_p") %in% methods)) {
    z_null <- .z_matrix(cumT_null, cumC_null, nn_looks)
    z_alt  <- .z_matrix(cumT_alt,  cumC_alt,  nn_looks)
  }

  # Calibrate OBF-style constant c if needed:
  # Reject if z_k >= c / sqrt(info_frac_k), i.e. max_k z_k * sqrt(info_frac_k) >= c.
  if ("gs_obf" %in% methods) {
    m_null <- apply(z_null * matrix(sqrt(info_frac), nrow = nrep, ncol = length(info_frac), byrow = TRUE),
                    1, max)
    gs_c <- as.numeric(stats::quantile(m_null, probs = 1 - alpha, names = FALSE))
    obf_bounds <- gs_c / sqrt(info_frac)
  }

  # --- Calibrate Bayesian threshold (if needed) ---
  if ("cal_bayes" %in% methods) {
    bayes_thresh <- .calibrate_bayes_threshold_vec(
      cumT_null, cumC_null, look_times, alpha,
      nrep_cal = min(nrep, 5000),
      posterior_method = posterior_method
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
        z <- if (scenario == "null") z_null else z_alt
        bounds_mat <- matrix(obf_bounds, nrow = nrep, ncol = length(look_times), byrow = TRUE)
        crossed <- z >= bounds_mat
        .store_results(crossed, look_times, Nmax, scenario, nrep,
                       rej_null, rej_alt, stop_null, stop_alt) -> tmp
        rej_null <- tmp$rej_null; rej_alt <- tmp$rej_alt
        stop_null <- tmp$stop_null; stop_alt <- tmp$stop_alt

      } else if (m == "naive_p") {
        z <- if (scenario == "null") z_null else z_alt
        pvals <- stats::pnorm(z, lower.tail = FALSE)
        crossed <- pvals < alpha
        .store_results(crossed, look_times, Nmax, scenario, nrep,
                       rej_null, rej_alt, stop_null, stop_alt) -> tmp
        rej_null <- tmp$rej_null; rej_alt <- tmp$rej_alt
        stop_null <- tmp$stop_null; stop_alt <- tmp$stop_alt

      } else if (m == "naive_bayes") {
        # Naive posterior threshold (not calibrated for optional stopping)
        pp <- .posterior_prob_matrix(cumT, cumC, nn_looks, method = posterior_method)
        crossed <- pp >= (1 - alpha)
        .store_results(crossed, look_times, Nmax, scenario, nrep,
                       rej_null, rej_alt, stop_null, stop_alt) -> tmp
        rej_null <- tmp$rej_null; rej_alt <- tmp$rej_alt
        stop_null <- tmp$stop_null; stop_alt <- tmp$stop_alt

      } else if (m == "cal_bayes") {
        # P(p_T > p_C | data) for independent Beta posteriors
        pp <- .posterior_prob_matrix(cumT, cumC, nn_looks, method = posterior_method)
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
                    lambda = lambda,
                    gs_c = gs_c,
                    bayes_thresh = if (exists("bayes_thresh")) bayes_thresh else NULL),
      raw = raw
    ),
    class = "ecomparison"
  )
}

# --- Internal helpers ---

#' Compute Wald z-statistics at each look (vectorized)
#' @keywords internal
.z_matrix <- function(cumT, cumC, nn) {
  phat_T <- cumT / nn
  phat_C <- cumC / nn
  delta <- phat_T - phat_C
  se <- sqrt(phat_T * (1 - phat_T) / nn + phat_C * (1 - phat_C) / nn)
  se[se < 1e-12] <- 1e-6
  delta / se
}

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

#' P(p_T > p_C | data) for Independent Beta Posteriors (Vectorized)
#'
#' Computes the posterior probability that the treatment arm response rate
#' exceeds the control arm response rate, given independent Beta posteriors
#' under a Jeffreys \code{Beta(0.5, 0.5)} prior.
#'
#' Two methods are available:
#' \describe{
#'   \item{\code{"normal"}}{Normal approximation to the Beta posterior
#'     difference. Each Beta is approximated by a Gaussian with matched
#'     mean and variance, so \eqn{P(p_T > p_C) \approx \Phi(\mu / \sigma)}
#'     where \eqn{\mu = E[p_T] - E[p_C]} and
#'     \eqn{\sigma = \sqrt{Var(p_T) + Var(p_C)}}. Extremely fast
#'     (fully vectorized, no loops) and accurate for per-arm \eqn{n \ge 10}.}
#'   \item{\code{"exact"}}{Exact computation via numerical integration of
#'     the regularized incomplete beta function:
#'     \deqn{P(p_T > p_C) = \int_0^1 I_x(a_C, b_C)\, f_{a_T,b_T}(x)\,dx}
#'     where \eqn{I_x(a, b) = } \code{pbeta(x, a, b)} is the regularized
#'     incomplete beta function and \eqn{f} is the Beta density. Exact to
#'     machine precision but requires element-wise \code{integrate} calls.}
#' }
#'
#' @param cumT Matrix (nrep x n_looks) of cumulative treatment successes.
#' @param cumC Matrix (nrep x n_looks) of cumulative control successes.
#' @param nn Matrix (nrep x n_looks) of per-arm sample sizes at each look.
#' @param method Character: \code{"normal"} (default) or \code{"exact"}.
#'
#' @return Matrix of same dimensions as \code{cumT} with posterior
#'   probabilities.
#'
#' @keywords internal
.posterior_prob_matrix <- function(cumT, cumC, nn, method = "normal") {
  # Jeffreys prior: Beta(0.5, 0.5)
  a_T <- 0.5 + cumT
  b_T <- 0.5 + nn - cumT
  a_C <- 0.5 + cumC
  b_C <- 0.5 + nn - cumC

  if (method == "exact") {
    # Exact computation using the regularized incomplete beta function.
    # P(p_T > p_C) = integral_0^1 pbeta(x, a_C, b_C) * dbeta(x, a_T, b_T) dx
    result <- matrix(NA_real_, nrow = nrow(cumT), ncol = ncol(cumT))
    for (j in seq_len(ncol(cumT))) {
      result[, j] <- mapply(
        function(aT, bT, aC, bC) {
          stats::integrate(
            function(x) stats::pbeta(x, aC, bC) * stats::dbeta(x, aT, bT),
            lower = 0, upper = 1,
            rel.tol = 1e-10
          )$value
        },
        a_T[, j], b_T[, j], a_C[, j], b_C[, j]
      )
    }
    return(result)
  }

  # Default: Normal approximation to Beta posteriors
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
#' @param posterior_method Passed to \code{.posterior_prob_matrix}.
#' @keywords internal
.calibrate_bayes_threshold_vec <- function(cumT_null, cumC_null, look_times,
                                           alpha, nrep_cal = 5000,
                                           posterior_method = "normal") {
  # Use at most nrep_cal rows from the null data
  nr <- min(nrow(cumT_null), nrep_cal)
  cumT <- cumT_null[seq_len(nr), , drop = FALSE]
  cumC <- cumC_null[seq_len(nr), , drop = FALSE]
  nn <- matrix(look_times, nrow = nr, ncol = length(look_times), byrow = TRUE)

  pp <- .posterior_prob_matrix(cumT, cumC, nn, method = posterior_method)
  max_pp <- apply(pp, 1, max)
  stats::quantile(max_pp, 1 - alpha, names = FALSE)
}

#' @keywords internal
.obf_boundaries <- function(look_times, Nmax, alpha) {
  # Kept for backward compatibility; the main simulation uses calibration
  info_frac <- look_times / Nmax
  z_alpha <- stats::qnorm(1 - alpha)
  z_alpha / sqrt(info_frac)
}

#' Scalar P(p_T > p_C) for Independent Beta Posteriors
#'
#' Computes \eqn{P(p_T > p_C)} for scalar Beta posterior parameters.
#' Supports both a fast normal approximation and the exact computation
#' via the regularized incomplete beta function.
#'
#' @param a_T,b_T Shape parameters of the treatment Beta posterior.
#' @param a_C,b_C Shape parameters of the control Beta posterior.
#' @param method Character: \code{"normal"} (default) for fast Gaussian
#'   approximation, or \code{"exact"} for exact numerical integration using
#'   the regularized incomplete beta function
#'   \eqn{I_x(a,b) = } \code{pbeta(x, a, b)}.
#'
#' @return Scalar probability.
#' @keywords internal
.posterior_prob_greater <- function(a_T, b_T, a_C, b_C,
                                   method = "normal") {
  if (method == "exact") {
    # Exact: P(p_T > p_C) = integral_0^1 pbeta(x, a_C, b_C) * dbeta(x, a_T, b_T) dx
    # Uses the regularized incomplete beta function I_x(a_C, b_C) = pbeta(x, a_C, b_C)
    return(stats::integrate(
      function(x) stats::pbeta(x, a_C, b_C) * stats::dbeta(x, a_T, b_T),
      lower = 0, upper = 1,
      rel.tol = 1e-10
    )$value)
  }

  # Normal approximation (default)
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
