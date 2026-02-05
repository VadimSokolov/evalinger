#' Design Calibration for E-Value Monitoring of Binary Trials
#'
#' Compute the optimal betting fraction, required maximum sample size, expected
#' stopping time, and power for a two-arm binary trial with e-value monitoring.
#'
#' @param p_C Control arm response probability.
#' @param delta Treatment effect (p_T - p_C) under the design alternative.
#' @param alpha Significance level (default 0.025).
#' @param power Target power (default 0.80).
#' @param lambda Optional fixed betting fraction. If NULL, uses the GROW-optimal value.
#' @param Nmax_max Upper bound for the sample size search (default 1000).
#' @param nrep Number of Monte Carlo replications for power estimation (default 10000).
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A list with components:
#' \describe{
#'   \item{p_C}{Control arm rate.}
#'   \item{p_T}{Treatment arm rate under the alternative.}
#'   \item{delta}{Treatment effect.}
#'   \item{lambda}{Betting fraction used.}
#'   \item{alpha}{Significance level.}
#'   \item{growth_rate}{Expected per-observation log-growth rate under the alternative.}
#'   \item{approx_stopping_time}{Approximate expected stopping time (analytic).}
#'   \item{Nmax_for_power}{Minimum Nmax per arm achieving the target power (by simulation).}
#'   \item{simulated_power}{Simulated power at \code{Nmax_for_power}.}
#'   \item{simulated_type1}{Simulated Type I error at \code{Nmax_for_power}.}
#' }
#'
#' @examples
#' des <- edesign_binary(p_C = 0.30, delta = 0.15, alpha = 0.025, power = 0.80,
#'                       nrep = 2000)
#' des
#'
#' @export
edesign_binary <- function(p_C, delta, alpha = 0.025, power = 0.80,
                           lambda = NULL, Nmax_max = 1000,
                           nrep = 10000, seed = 42) {
  p_T <- p_C + delta
  stopifnot(p_T > p_C, p_T < 1, p_C > 0, p_C < 1)
  stopifnot(alpha > 0, alpha < 1, power > 0, power < 1)

  if (is.null(lambda)) {
    lambda <- grow_lambda(p_T, p_C)
  }

  g <- expected_growth_rate(lambda, p_T, p_C)
  approx_tau <- if (g > 0) log(1 / alpha) / g else Inf

  # Binary search for Nmax that achieves target power
  threshold <- log(1 / alpha)
  set.seed(seed)

  sim_power <- function(Nmax) {
    rng <- .Random.seed
    rejected <- 0L
    for (i in seq_len(nrep)) {
      x_T <- stats::rbinom(Nmax, 1, p_T)
      x_C <- stats::rbinom(Nmax, 1, p_C)
      D <- x_T - x_C
      log_e <- cumsum(log(pmax(1 + lambda * D, .Machine$double.xmin)))
      if (any(log_e >= threshold)) rejected <- rejected + 1L
    }
    rejected / nrep
  }

  sim_type1 <- function(Nmax) {
    rejected <- 0L
    for (i in seq_len(nrep)) {
      x_T <- stats::rbinom(Nmax, 1, p_C)  # null: both at p_C
      x_C <- stats::rbinom(Nmax, 1, p_C)
      D <- x_T - x_C
      log_e <- cumsum(log(pmax(1 + lambda * D, .Machine$double.xmin)))
      if (any(log_e >= threshold)) rejected <- rejected + 1L
    }
    rejected / nrep
  }

  # Search: find smallest Nmax in multiples of 10 that achieves power
  Nmax_found <- Nmax_max
  for (Nmax_candidate in seq(20, Nmax_max, by = 10)) {
    set.seed(seed)
    pw <- sim_power(Nmax_candidate)
    if (pw >= power) {
      Nmax_found <- Nmax_candidate
      break
    }
  }

  set.seed(seed)
  final_power <- sim_power(Nmax_found)
  set.seed(seed + 1)
  final_type1 <- sim_type1(Nmax_found)

  structure(
    list(
      p_C = p_C, p_T = p_T, delta = delta,
      lambda = lambda, alpha = alpha,
      growth_rate = g,
      approx_stopping_time = approx_tau,
      Nmax_for_power = Nmax_found,
      simulated_power = final_power,
      simulated_type1 = final_type1
    ),
    class = "edesign"
  )
}

#' Print an E-Value Design Object
#'
#' Display the design parameters for an e-value-based trial, including
#' arm response rates, betting fraction, growth rate, expected stopping
#' time, and simulated operating characteristics.
#'
#' @param x An \code{"edesign"} object from \code{\link{edesign_binary}}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.edesign <- function(x, ...) {
  cat("E-value design for two-arm binary trial\n")
  cat("========================================\n")
  cat("  Control rate (p_C):     ", x$p_C, "\n")
  cat("  Treatment rate (p_T):   ", x$p_T, "\n")
  cat("  Effect size (delta):    ", x$delta, "\n")
  cat("  Betting fraction:       ", sprintf("%.4f", x$lambda), "\n")
  cat("  Alpha:                  ", x$alpha, "\n")
  cat("  Growth rate g(lambda):  ", sprintf("%.6f", x$growth_rate), "\n")
  cat("  Approx. stopping time:  ", sprintf("%.0f", x$approx_stopping_time), " per arm\n")
  cat("  Nmax for target power:  ", x$Nmax_for_power, " per arm\n")
  cat("  Simulated power:        ", sprintf("%.3f", x$simulated_power), "\n")
  cat("  Simulated Type I error: ", sprintf("%.4f", x$simulated_type1), "\n")
  invisible(x)
}
