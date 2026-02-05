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
#' posterior probability \eqn{P(p_T > p_C \mid \text{data})}, and calibrates
#' the decision threshold to achieve the nominal Type I error rate.
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

  # OBF spending function: alpha_k = 2 - 2*Phi(z_alpha / sqrt(t_k))
  obf_bounds <- .obf_boundaries(look_times, Nmax, alpha)

  # Calibrated Bayes: find posterior threshold for Type I error control
  bayes_thresh <- .calibrate_bayes_threshold(p_C, look_times, alpha,
                                             nrep_cal = min(nrep, 5000),
                                             seed = seed + 999)

  # Pre-allocate results
  raw <- list()
  for (m in methods) {
    raw[[m]] <- list(
      null_rej = logical(nrep), alt_rej = logical(nrep),
      null_stop = integer(nrep), alt_stop = integer(nrep)
    )
  }

  set.seed(seed)
  for (rep in seq_len(nrep)) {
    # Generate full data (both null and alt)
    x_C_null <- stats::rbinom(Nmax, 1, p_C)
    x_T_null <- stats::rbinom(Nmax, 1, p_C)
    x_C_alt  <- stats::rbinom(Nmax, 1, p_C)
    x_T_alt  <- stats::rbinom(Nmax, 1, p_T_alt)

    for (m in methods) {
      for (scenario in c("null", "alt")) {
        x_T <- if (scenario == "null") x_T_null else x_T_alt
        x_C <- if (scenario == "null") x_C_null else x_C_alt

        result <- .run_method(m, x_T, x_C, look_times, lambda,
                              threshold_log_e, obf_bounds, bayes_thresh,
                              alpha, Nmax)

        tag_rej <- paste0(scenario, "_rej")
        tag_stop <- paste0(scenario, "_stop")
        raw[[m]][[tag_rej]][rep] <- result$rejected
        raw[[m]][[tag_stop]][rep] <- result$stop_n
      }
    }
  }

  # Summarize
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
    if (log_e[t] >= threshold) {
      return(list(rejected = TRUE, stop_n = t))
    }
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.run_gs <- function(x_T, x_C, look_times, obf_bounds, Nmax) {
  for (k in seq_along(look_times)) {
    t <- look_times[k]
    n_T <- sum(x_T[1:t])
    n_C <- sum(x_C[1:t])
    p_hat_T <- n_T / t
    p_hat_C <- n_C / t
    delta <- p_hat_T - p_hat_C
    se <- sqrt(p_hat_T * (1 - p_hat_T) / t + p_hat_C * (1 - p_hat_C) / t)
    if (se < 1e-12) se <- 1e-6
    z <- delta / se
    if (z >= obf_bounds[k]) {
      return(list(rejected = TRUE, stop_n = t))
    }
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.run_naive_p <- function(x_T, x_C, look_times, alpha, Nmax) {
  for (t in look_times) {
    n_T <- sum(x_T[1:t])
    n_C <- sum(x_C[1:t])
    p_hat_T <- n_T / t
    p_hat_C <- n_C / t
    delta <- p_hat_T - p_hat_C
    se <- sqrt(p_hat_T * (1 - p_hat_T) / t + p_hat_C * (1 - p_hat_C) / t)
    if (se < 1e-12) se <- 1e-6
    z <- delta / se
    pval <- stats::pnorm(z, lower.tail = FALSE)
    if (pval < alpha) {
      return(list(rejected = TRUE, stop_n = t))
    }
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.run_bayes <- function(x_T, x_C, look_times, threshold, Nmax) {
  for (t in look_times) {
    n_T <- sum(x_T[1:t])
    n_C <- sum(x_C[1:t])
    # Jeffreys prior: Beta(0.5, 0.5)
    a_T <- 0.5 + n_T
    b_T <- 0.5 + t - n_T
    a_C <- 0.5 + n_C
    b_C <- 0.5 + t - n_C
    # P(p_T > p_C | data) by Monte Carlo
    pp <- .posterior_prob_greater(a_T, b_T, a_C, b_C, nsim = 5000)
    if (pp >= threshold) {
      return(list(rejected = TRUE, stop_n = t))
    }
  }
  list(rejected = FALSE, stop_n = Nmax)
}

#' @keywords internal
.posterior_prob_greater <- function(a_T, b_T, a_C, b_C, nsim = 5000) {
  samp_T <- stats::rbeta(nsim, a_T, b_T)
  samp_C <- stats::rbeta(nsim, a_C, b_C)
  mean(samp_T > samp_C)
}

#' @keywords internal
.obf_boundaries <- function(look_times, Nmax, alpha) {
  # O'Brien-Fleming-like spending: boundary at look k is z_alpha / sqrt(t_k)
  # where t_k = look_times[k] / Nmax (information fraction)
  info_frac <- look_times / Nmax
  z_alpha <- stats::qnorm(1 - alpha)
  z_alpha / sqrt(info_frac)
}

#' @keywords internal
.calibrate_bayes_threshold <- function(p_C, look_times, alpha,
                                       nrep_cal = 5000, seed = 999) {
  # Find the posterior probability threshold that controls Type I error at alpha
  # under sequential monitoring with Jeffreys prior
  set.seed(seed)
  Nmax <- max(look_times)

  # Simulate null trials
  max_pp <- numeric(nrep_cal)
  for (rep in seq_len(nrep_cal)) {
    x_T <- stats::rbinom(Nmax, 1, p_C)
    x_C <- stats::rbinom(Nmax, 1, p_C)
    pp_vec <- numeric(length(look_times))
    for (k in seq_along(look_times)) {
      t <- look_times[k]
      n_T <- sum(x_T[1:t])
      n_C <- sum(x_C[1:t])
      a_T <- 0.5 + n_T
      b_T <- 0.5 + t - n_T
      a_C <- 0.5 + n_C
      b_C <- 0.5 + t - n_C
      pp_vec[k] <- .posterior_prob_greater(a_T, b_T, a_C, b_C, nsim = 2000)
    }
    max_pp[rep] <- max(pp_vec)
  }

  # Threshold: 1 - alpha quantile of the max posterior probability under null
  thresh <- stats::quantile(max_pp, 1 - alpha, names = FALSE)
  thresh
}
