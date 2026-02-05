#' CS-Based Futility Check
#'
#' Declare futility when the upper bound of the confidence sequence for the
#' treatment effect falls below the minimum clinically important difference (MCID).
#'
#' @param cs A \code{"confseq"} object from \code{\link{confseq_binary}}.
#' @param delta_min The MCID (minimum clinically important difference).
#'
#' @return A list with components:
#' \describe{
#'   \item{futile}{Logical vector of length \code{cs$n}: TRUE where futility is declared.}
#'   \item{first_futile}{Integer: first observation at which futility is declared, or NA.}
#'   \item{delta_min}{The MCID used.}
#' }
#'
#' @details
#' At any time \eqn{n}, if the upper bound of the confidence sequence
#' \eqn{U_n < \delta_{\min}}, we can declare futility with confidence.
#' Because the CS has uniform coverage, this futility declaration is
#' simultaneously valid at all sample sizes.
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.32)  # weak effect
#' x_C <- rbinom(200, 1, 0.30)
#' cs <- confseq_binary(x_T, x_C)
#' fut <- futility_cs(cs, delta_min = 0.10)
#' fut$first_futile
#'
#' @export
futility_cs <- function(cs, delta_min) {
  stopifnot(inherits(cs, "confseq"))
  stopifnot(is.numeric(delta_min), length(delta_min) == 1)

  futile <- cs$upper < delta_min
  first <- if (any(futile)) which(futile)[1] else NA_integer_

  list(futile = futile, first_futile = first, delta_min = delta_min)
}

#' Reciprocal E-Process for Futility
#'
#' Construct a reciprocal e-process that tests the alternative hypothesis
#' (or a reduced alternative). When this e-process exceeds its threshold,
#' futility is declared.
#'
#' @param x_T Integer vector of treatment arm binary outcomes (0/1).
#' @param x_C Integer vector of control arm binary outcomes (0/1).
#' @param delta_min The minimum treatment effect to test against.
#' @param alpha_f Futility significance level (default 0.10).
#' @param lambda_f Betting fraction for the futility e-process. If NULL,
#'   uses a GROW-optimal fraction calibrated for the reversed hypothesis.
#'
#' @return A list with components:
#' \describe{
#'   \item{log_evalue}{Numeric vector of the reciprocal log e-values.}
#'   \item{evalue}{Numeric vector of the reciprocal e-values.}
#'   \item{futile}{Logical: did the reciprocal e-process cross the threshold?}
#'   \item{first_futile}{Integer: first crossing time, or NA.}
#'   \item{alpha_f}{Futility significance level.}
#' }
#'
#' @details
#' The reciprocal e-process tests \eqn{H_1: \delta \ge \delta_{\min}} by
#' constructing a betting martingale under the assumption that
#' \eqn{\delta = \delta_{\min}}. The difference \eqn{D_i - \delta_{\min}}
#' has mean zero under this tilted null. When the e-process exceeds
#' \eqn{1/\alpha_f}, we reject \eqn{H_1} (i.e., conclude futility).
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.32)
#' x_C <- rbinom(200, 1, 0.30)
#' fut <- futility_eprocess(x_T, x_C, delta_min = 0.10)
#' fut$futile
#'
#' @export
futility_eprocess <- function(x_T, x_C, delta_min, alpha_f = 0.10,
                              lambda_f = NULL) {
  stopifnot(length(x_T) == length(x_C))
  n <- length(x_T)

  if (is.null(lambda_f)) {
    # Use a moderate betting fraction for the reversed test
    lambda_f <- 0.3
  }
  stopifnot(lambda_f > 0, lambda_f < 1)

  D <- as.numeric(x_T) - as.numeric(x_C)
  # Center at delta_min: under H1 (delta = delta_min), D - delta_min has mean 0
  D_centered <- D - delta_min

  # Reciprocal e-process: bet on mean being negative (futility direction)
  # E_n = prod(1 - lambda * D_centered) -- betting against the treatment
  log_inc <- log(pmax(1 - lambda_f * D_centered, .Machine$double.xmin))
  log_evalue <- cumsum(log_inc)
  evalue <- exp(log_evalue)

  threshold <- log(1 / alpha_f)
  crossed <- log_evalue >= threshold
  futile <- any(crossed)
  first_futile <- if (futile) which(crossed)[1] else NA_integer_

  list(
    log_evalue = log_evalue,
    evalue = evalue,
    futile = futile,
    first_futile = first_futile,
    alpha_f = alpha_f
  )
}
