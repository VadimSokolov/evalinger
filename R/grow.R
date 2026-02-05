#' GROW Optimal Betting Fraction
#'
#' Compute the Growth Rate Optimal in the Worst case (GROW) betting fraction
#' for a two-arm binary trial. This is the value of lambda that maximizes
#' the expected log-growth rate of the betting martingale under a specified
#' alternative.
#'
#' @param p_T Response probability under the alternative for the treatment arm.
#' @param p_C Response probability for the control arm (same under null and alternative).
#'
#' @return A scalar betting fraction in (0, 1).
#'
#' @details
#' For the betting martingale \eqn{E_n = \prod_{i=1}^n (1 + \lambda D_i)},
#' where \eqn{D_i = X_i^T - X_i^C \in \{-1, 0, 1\}}, the expected
#' per-observation log-growth rate under \eqn{(p_T, p_C)} is
#'
#' \deqn{g(\lambda) = p_T(1-p_C)\log(1+\lambda) + (1-p_T)p_C\log(1-\lambda).}
#'
#' The GROW-optimal fraction is
#'
#' \deqn{\lambda^* = \frac{p_T(1-p_C) - (1-p_T)p_C}{p_T(1-p_C) + (1-p_T)p_C}.}
#'
#' @examples
#' # Design alternative: p_T = 0.45, p_C = 0.30
#' grow_lambda(0.45, 0.30)
#'
#' @export
grow_lambda <- function(p_T, p_C) {
  stopifnot(
    is.numeric(p_T), length(p_T) == 1, p_T > 0, p_T < 1,
    is.numeric(p_C), length(p_C) == 1, p_C > 0, p_C < 1,
    p_T > p_C
  )
  a <- p_T * (1 - p_C)
  b <- (1 - p_T) * p_C
  lam <- (a - b) / (a + b)
  lam
}

#' Expected Log-Growth Rate of the Betting Martingale
#'
#' Compute the expected per-observation log-growth rate \eqn{g(\lambda)}
#' of the betting martingale under a specified alternative.
#'
#' @param lambda Betting fraction in (0, 1).
#' @param p_T Treatment arm response probability.
#' @param p_C Control arm response probability.
#'
#' @return A scalar: the expected log-growth rate per observation.
#'
#' @examples
#' lam <- grow_lambda(0.45, 0.30)
#' expected_growth_rate(lam, 0.45, 0.30)
#'
#' @export
expected_growth_rate <- function(lambda, p_T, p_C) {
  stopifnot(
    is.numeric(lambda), length(lambda) == 1,
    lambda > 0, lambda < 1,
    is.numeric(p_T), is.numeric(p_C),
    p_T > 0, p_T < 1, p_C > 0, p_C < 1
  )
  a <- p_T * (1 - p_C)
  b <- (1 - p_T) * p_C
  g <- a * log(1 + lambda) + b * log(1 - lambda)
  g
}

#' Expected Stopping Time of the Betting Martingale
#'
#' Approximate the expected number of observations (per arm) for the
#' e-process to cross the rejection threshold \eqn{1/\alpha}.
#'
#' @param lambda Betting fraction in (0, 1).
#' @param p_T Treatment arm response probability under the alternative.
#' @param p_C Control arm response probability.
#' @param alpha Significance level (default 0.025).
#'
#' @return A scalar: the approximate expected stopping time (per arm).
#'
#' @details
#' Uses the approximation \eqn{\tau \approx \log(1/\alpha) / g(\lambda)},
#' where \eqn{g} is the expected log-growth rate. This is a first-order
#' approximation that tends to underestimate the true expected stopping time.
#'
#' @examples
#' lam <- grow_lambda(0.45, 0.30)
#' expected_stopping_time(lam, 0.45, 0.30, alpha = 0.025)
#'
#' @export
expected_stopping_time <- function(lambda, p_T, p_C, alpha = 0.025) {
  g <- expected_growth_rate(lambda, p_T, p_C)
  if (g <= 0) return(Inf)
  log(1 / alpha) / g
}

#' Grid Search over Betting Fractions
#'
#' Evaluate the expected growth rate and stopping time for a grid of
#' lambda values and design alternatives, useful for sensitivity analysis.
#'
#' @param p_C Control arm response probability.
#' @param delta_grid Numeric vector of treatment effects (p_T - p_C) to evaluate.
#' @param lambda_grid Numeric vector of betting fractions to evaluate.
#'   Defaults to \code{seq(0.05, 0.95, by = 0.05)}.
#' @param alpha Significance level (default 0.025).
#'
#' @return A data frame with columns \code{delta}, \code{lambda},
#'   \code{growth_rate}, and \code{expected_n}.
#'
#' @examples
#' grow_lambda_grid(p_C = 0.30, delta_grid = c(0.10, 0.15, 0.20))
#'
#' @export
grow_lambda_grid <- function(p_C, delta_grid,
                             lambda_grid = seq(0.05, 0.95, by = 0.05),
                             alpha = 0.025) {
  stopifnot(is.numeric(p_C), length(p_C) == 1, p_C > 0, p_C < 1)
  stopifnot(is.numeric(delta_grid), all(delta_grid > 0))
  stopifnot(is.numeric(lambda_grid), all(lambda_grid > 0), all(lambda_grid < 1))

  rows <- vector("list", length(delta_grid) * length(lambda_grid))
  k <- 1L
  for (delta in delta_grid) {
    p_T <- p_C + delta
    if (p_T >= 1) next
    for (lam in lambda_grid) {
      g <- expected_growth_rate(lam, p_T, p_C)
      en <- if (g > 0) log(1 / alpha) / g else Inf
      rows[[k]] <- data.frame(
        delta = delta, lambda = lam,
        growth_rate = g, expected_n = en,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, rows[seq_len(k - 1L)])
}
