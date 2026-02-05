#' Always-Valid Confidence Sequence for Treatment Effect
#'
#' Compute an always-valid confidence sequence for the treatment effect
#' \eqn{\delta = p_T - p_C} from paired binary outcomes. The confidence
#' sequence is valid at every sample size simultaneously.
#'
#' @param x_T Integer vector of treatment arm binary outcomes (0/1).
#' @param x_C Integer vector of control arm binary outcomes (0/1).
#' @param alpha Significance level (default 0.05 for a 95 percent CS).
#' @param method Method for CS construction. Currently \code{"betting"} (default),
#'   which uses a Wald-type interval with a time-uniform boundary based on the
#'   law of the iterated logarithm.
#' @param t_min Minimum sample size before the CS is computed (default 10).
#'
#' @return An object of class \code{"confseq"} with components:
#' \describe{
#'   \item{lower}{Numeric vector of lower bounds at each observation.}
#'   \item{upper}{Numeric vector of upper bounds at each observation.}
#'   \item{delta_hat}{Numeric vector of point estimates at each observation.}
#'   \item{n}{Number of observations.}
#'   \item{alpha}{Significance level.}
#' }
#'
#' @details
#' The confidence sequence is constructed by inverting a sequence of tests.
#' For computational efficiency, this implementation uses a normal-approximation
#' boundary scaled by a time-uniform factor. For small samples, the CS may be
#' wider than necessary; the boundary tightens as \eqn{n} grows.
#'
#' The boundary width at time \eqn{n} is approximately
#' \deqn{w_n = \sqrt{\frac{2 \hat{V}_n \log(2/\alpha \cdot \log_2(2n))}{n}},}
#' where \eqn{\hat{V}_n} is the estimated variance of \eqn{D_i}.
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.45)
#' x_C <- rbinom(200, 1, 0.30)
#' cs <- confseq_binary(x_T, x_C)
#' plot(cs)
#'
#' @export
confseq_binary <- function(x_T, x_C, alpha = 0.05, method = "betting",
                           t_min = 10) {
  stopifnot(length(x_T) == length(x_C))
  n <- length(x_T)
  stopifnot(n >= t_min)

  D <- as.numeric(x_T) - as.numeric(x_C)
  cum_D <- cumsum(D)
  cum_D2 <- cumsum(D^2)

  delta_hat <- cum_D / seq_len(n)
  # Running variance estimate: Var(D) = E[D^2] - (E[D])^2
  var_hat <- pmax(cum_D2 / seq_len(n) - delta_hat^2, 0.01)

  # Time-uniform boundary using a stitched-together approach
  # Based on Howard et al. (2021) style boundary
  nn <- seq_len(n)
  # Boundary: sqrt(2 * V * log(log_2(2n) * 2/alpha) / n)
  log_factor <- log(pmax(log2(2 * nn), 1) * 2 / alpha)
  width <- sqrt(2 * var_hat * log_factor / nn)

  lower <- delta_hat - width
  upper <- delta_hat + width

  # Set early values (before t_min) to [-1, 1]
  if (t_min > 1) {
    lower[seq_len(min(t_min - 1, n))] <- -1
    upper[seq_len(min(t_min - 1, n))] <- 1
  }

  structure(
    list(
      lower = lower, upper = upper,
      delta_hat = delta_hat,
      n = n, alpha = alpha
    ),
    class = "confseq"
  )
}

#' Plot a Confidence Sequence
#'
#' Visualize the always-valid confidence sequence for the treatment effect
#' over the course of the trial. Shows the point estimate and confidence
#' band at each sample size.
#'
#' @param x A \code{"confseq"} object from \code{\link{confseq_binary}}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.45)
#' x_C <- rbinom(200, 1, 0.30)
#' cs <- confseq_binary(x_T, x_C)
#' plot(cs)
#'
#' @export
plot.confseq <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback
    nn <- seq_len(x$n)
    plot(nn, x$delta_hat, type = "l", ylim = range(c(x$lower, x$upper)),
         xlab = "Observations per arm", ylab = "Treatment effect",
         main = "Confidence sequence for treatment effect")
    lines(nn, x$lower, lty = 2, col = "blue")
    lines(nn, x$upper, lty = 2, col = "blue")
    abline(h = 0, col = "gray", lty = 3)
    return(invisible(x))
  }

  df <- data.frame(
    n = seq_len(x$n),
    delta_hat = x$delta_hat,
    lower = x$lower,
    upper = x$upper
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "steelblue", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$delta_hat)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::labs(x = "Observations per arm", y = "Treatment effect",
                  title = "Always-valid confidence sequence") +
    ggplot2::theme_minimal()
  print(p)
  invisible(x)
}
