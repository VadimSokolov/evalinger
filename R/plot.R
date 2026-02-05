#' Plot an E-Process
#'
#' Visualize the sample path of an e-process with the rejection threshold.
#'
#' @param x An \code{"eprocess"} object.
#' @param log_scale Logical: plot on the log scale? Default TRUE.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the \code{ggplot} object if \pkg{ggplot2} is
#'   available, or the base-R plot otherwise.
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.45)
#' x_C <- rbinom(200, 1, 0.30)
#' ep <- eprocess_binary(x_T, x_C, lambda = 0.31)
#' plot(ep)
#'
#' @export
plot.eprocess <- function(x, log_scale = TRUE, ...) {
  nn <- seq_len(x$n)
  y <- if (log_scale) x$log_evalue else x$evalue
  ylab <- if (log_scale) "log(E)" else "E-value"
  thresh <- if (log_scale) x$threshold else exp(x$threshold)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot(nn, y, type = "l", col = "steelblue", lwd = 1.5,
         xlab = "Observations per arm", ylab = ylab,
         main = "E-process sample path")
    abline(h = thresh, col = "red", lty = 2, lwd = 1.5)
    if (x$rejected) {
      abline(v = x$rejection_time, col = "gray50", lty = 3)
    }
    legend("topleft",
           legend = c("E-process", paste0("Threshold (1/alpha = ",
                      sprintf("%.0f", 1/x$alpha), ")")),
           col = c("steelblue", "red"), lty = c(1, 2), lwd = 1.5,
           bty = "n")
    return(invisible(x))
  }

  df <- data.frame(n = nn, y = y)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$n, y = .data$y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = thresh, linetype = "dashed",
                        color = "red", linewidth = 0.7) +
    ggplot2::labs(x = "Observations per arm", y = ylab,
                  title = "E-process sample path") +
    ggplot2::theme_minimal(base_size = 12)

  if (x$rejected) {
    p <- p + ggplot2::geom_vline(xintercept = x$rejection_time,
                                 linetype = "dotted", color = "gray50")
  }

  print(p)
  invisible(p)
}

#' Plot Hybrid E-Process + Group Sequential Overlay
#'
#' Produce a dual-panel or overlaid visualization of an e-process alongside
#' a group sequential boundary, for comparing the two monitoring approaches
#' on the same data.
#'
#' @param eproc An \code{"eprocess"} object.
#' @param look_times Integer vector of interim analysis times.
#' @param Nmax Maximum sample size per arm.
#' @param alpha Significance level (default 0.025).
#' @param x_T Integer vector of treatment arm outcomes.
#' @param x_C Integer vector of control arm outcomes.
#' @param gs_c Optional calibrated OBF constant \eqn{c} for an OBF-style boundary
#'   \eqn{z \ge c/\sqrt{t}}. If NULL (default), uses \code{qnorm(1 - alpha)}.
#'
#' @return Invisibly returns the \code{ggplot} object (or invisible NULL
#'   for base graphics).
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.45)
#' x_C <- rbinom(200, 1, 0.30)
#' ep <- eprocess_binary(x_T, x_C, lambda = 0.31)
#' plot_hybrid(ep, look_times = c(50, 100, 150, 200), Nmax = 200,
#'             x_T = x_T, x_C = x_C)
#'
#' @export
plot_hybrid <- function(eproc, look_times, Nmax, alpha = 0.025,
                        x_T, x_C, gs_c = NULL) {
  stopifnot(inherits(eproc, "eprocess"))
  info_frac <- look_times / Nmax
  z_alpha <- if (is.null(gs_c)) stats::qnorm(1 - alpha) else gs_c
  obf <- z_alpha / sqrt(info_frac)

  # Compute z-statistics at each look
  z_stats <- numeric(length(look_times))
  for (k in seq_along(look_times)) {
    t <- look_times[k]
    p_T <- sum(x_T[1:t]) / t
    p_C <- sum(x_C[1:t]) / t
    se <- sqrt(p_T * (1 - p_T) / t + p_C * (1 - p_C) / t)
    if (se < 1e-12) se <- 1e-6
    z_stats[k] <- (p_T - p_C) / se
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    # E-process panel
    nn <- seq_len(eproc$n)
    plot(nn, eproc$log_evalue, type = "l", col = "steelblue", lwd = 1.5,
         xlab = "n per arm", ylab = "log(E)", main = "E-process monitoring")
    abline(h = eproc$threshold, col = "red", lty = 2, lwd = 1.5)
    legend("topleft",
           legend = c("E-process", bquote("E-value threshold log(1/" * alpha * ")")),
           col = c("steelblue", "red"), lty = c(1, 2), lwd = 1.5, bty = "n")
    # GS panel
    plot(look_times, z_stats, type = "b", pch = 19, col = "darkorange",
         ylim = range(c(z_stats, obf)),
         xlab = "n per arm", ylab = "Z-statistic",
         main = "Group sequential (OBF)")
    lines(look_times, obf, col = "red", lty = 2, lwd = 1.5)
    legend("topright",
           legend = c("Z-statistic", "GS boundary (OBF)"),
           col = c("darkorange", "red"), lty = c(1, 2), pch = c(19, NA),
           lwd = 1.5, bty = "n")
    par(mfrow = c(1, 1))
    return(invisible(NULL))
  }

  df_e <- data.frame(n = seq_len(eproc$n), log_e = eproc$log_evalue)
  df_gs <- data.frame(n = look_times, z = z_stats, bound = obf)

  thresh_label <- paste0("E-value threshold log(1/a) = ",
                         sprintf("%.1f", eproc$threshold))

  p1 <- ggplot2::ggplot(df_e, ggplot2::aes(x = .data$n, y = .data$log_e)) +
    ggplot2::geom_line(ggplot2::aes(color = "E-process"), linewidth = 0.7) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = eproc$threshold,
                                     linetype = thresh_label),
                        color = "red", linewidth = 0.7) +
    ggplot2::scale_color_manual(name = NULL,
                                values = c("E-process" = "steelblue")) +
    ggplot2::scale_linetype_manual(name = NULL,
                                   values = stats::setNames("dashed", thresh_label)) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1,
        override.aes = list(linetype = "solid")),
      linetype = ggplot2::guide_legend(order = 2,
        override.aes = list(color = "red"))
    ) +
    ggplot2::labs(x = "Observations per arm", y = "log(E)",
                  title = "E-process monitoring") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  p2 <- ggplot2::ggplot(df_gs) +
    ggplot2::geom_point(ggplot2::aes(x = .data$n, y = .data$z,
                                     color = "Z-statistic"),
                        size = 3) +
    ggplot2::geom_line(ggplot2::aes(x = .data$n, y = .data$z,
                                    color = "Z-statistic")) +
    ggplot2::geom_line(ggplot2::aes(x = .data$n, y = .data$bound,
                                    color = "GS boundary (OBF)"),
                       linetype = "dashed") +
    ggplot2::scale_color_manual(name = NULL,
                                values = c("Z-statistic" = "darkorange",
                                           "GS boundary (OBF)" = "red")) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(linetype = c("dashed", "solid"),
                            shape = c(NA, 19)))
    ) +
    ggplot2::labs(x = "Observations per arm", y = "Z-statistic",
                  title = "Group sequential (OBF boundary)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")

  # Stack vertically if patchwork is available, otherwise print side by side
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p1 / p2
    print(combined)
    return(invisible(combined))
  }

  print(p1)
  print(p2)
  invisible(list(p1 = p1, p2 = p2))
}

#' Plot Comparison Results
#'
#' Bar chart showing Type I error, power, and average sample size across
#' monitoring methods.
#'
#' @param x An \code{"ecomparison"} object from \code{\link{simulate_comparison}}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the \code{ggplot} object or invisible \code{NULL}.
#'
#' @examples
#' cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
#'                            n_looks = 4, nrep = 500, seed = 1)
#' plot(cmp)
#'
#' @export
plot.ecomparison <- function(x, ...) {
  res <- x$results

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    par(mfrow = c(1, 3), mar = c(6, 4, 3, 1))
    barplot(res$null_rej, names.arg = res$method, main = "Type I Error",
            col = "salmon", las = 2, ylab = "Rejection rate")
    abline(h = x$design$alpha, col = "red", lty = 2)
    barplot(res$alt_rej, names.arg = res$method, main = "Power",
            col = "steelblue", las = 2, ylab = "Rejection rate")
    barplot(res$avg_n_alt, names.arg = res$method, main = "Avg N (alt)",
            col = "seagreen", las = 2, ylab = "Average N per arm")
    par(mfrow = c(1, 1))
    return(invisible(NULL))
  }

  # Type I error panel
  p1 <- ggplot2::ggplot(res, ggplot2::aes(x = .data$method, y = .data$null_rej)) +
    ggplot2::geom_col(fill = "salmon", width = 0.6) +
    ggplot2::geom_hline(yintercept = x$design$alpha, linetype = "dashed",
                        color = "red") +
    ggplot2::labs(x = NULL, y = "Rejection rate", title = "Type I error") +
    ggplot2::theme_minimal(base_size = 11)

  # Power panel
  p2 <- ggplot2::ggplot(res, ggplot2::aes(x = .data$method, y = .data$alt_rej)) +
    ggplot2::geom_col(fill = "steelblue", width = 0.6) +
    ggplot2::labs(x = NULL, y = "Rejection rate", title = "Power") +
    ggplot2::theme_minimal(base_size = 11)

  # Avg N panel
  p3 <- ggplot2::ggplot(res, ggplot2::aes(x = .data$method, y = .data$avg_n_alt)) +
    ggplot2::geom_col(fill = "seagreen", width = 0.6) +
    ggplot2::labs(x = NULL, y = "Average N per arm",
                  title = "Expected sample size (alternative)") +
    ggplot2::theme_minimal(base_size = 11)

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p1 + p2 + p3
    print(combined)
    return(invisible(combined))
  }

  print(p1)
  print(p2)
  print(p3)
  invisible(list(p1 = p1, p2 = p2, p3 = p3))
}

#' Plot Comparison Results (Convenience Wrapper)
#'
#' A convenience wrapper for \code{\link{plot.ecomparison}} that dispatches
#' to the S3 plot method for \code{"ecomparison"} objects.
#'
#' @param x An \code{"ecomparison"} object from \code{\link{simulate_comparison}}.
#' @param ... Additional arguments passed to \code{\link{plot.ecomparison}}.
#'
#' @return Invisibly returns the \code{ggplot} object(s) or invisible \code{NULL}.
#'
#' @examples
#' cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
#'                            n_looks = 4, nrep = 500, seed = 1)
#' plot_comparison(cmp)
#'
#' @seealso \code{\link{simulate_comparison}}, \code{\link{plot.ecomparison}}
#'
#' @export
plot_comparison <- function(x, ...) {
  if (!inherits(x, "ecomparison")) {
    stop("x must be an 'ecomparison' object from simulate_comparison()")
  }
  plot.ecomparison(x, ...)
}
