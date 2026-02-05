#' Print an E-Process Object
#'
#' Display a concise summary of an e-process, including the endpoint type,
#' sample size, betting fraction, final e-value, and rejection status.
#'
#' @param x An \code{"eprocess"} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.eprocess <- function(x, ...) {
  cat("E-process (", x$endpoint, " endpoint)\n", sep = "")
  cat("  Observations (per arm): ", x$n, "\n", sep = "")
  lam_str <- if (length(x$lambda) == 1) {
    sprintf("%.4f (fixed)", x$lambda)
  } else {
    sprintf("adaptive (%d values)", length(x$lambda))
  }
  cat("  Betting fraction: ", lam_str, "\n", sep = "")
  cat("  Final log(E):     ", sprintf("%.3f", x$log_evalue[x$n]), "\n", sep = "")
  cat("  Final E:          ", sprintf("%.3f", x$evalue[x$n]), "\n", sep = "")
  cat("  Threshold log(1/alpha): ", sprintf("%.3f", x$threshold), "\n", sep = "")
  if (x$rejected) {
    cat("  REJECTED at observation ", x$rejection_time, "\n", sep = "")
  } else {
    cat("  Not rejected\n")
  }
  invisible(x)
}

#' Summarize an E-Process Object
#'
#' Display a detailed summary of an e-process, including the endpoint type,
#' sample size, betting fraction, alpha level, final and maximum e-values,
#' always-valid p-value, and rejection status.
#'
#' @param object An \code{"eprocess"} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
summary.eprocess <- function(object, ...) {
  cat("E-process summary\n")
  cat("=================\n")
  cat("Endpoint type:      ", object$endpoint, "\n")
  cat("Observations:       ", object$n, "\n")
  lam_str <- if (length(object$lambda) == 1) {
    sprintf("%.4f", object$lambda)
  } else {
    sprintf("adaptive (mean %.4f)", mean(object$lambda))
  }
  cat("Betting fraction:   ", lam_str, "\n")
  cat("Alpha:              ", object$alpha, "\n")
  cat("Threshold log(1/a): ", sprintf("%.4f", object$threshold), "\n")
  cat("Final log(E):       ", sprintf("%.4f", object$log_evalue[object$n]), "\n")
  cat("Final E-value:      ", sprintf("%.4f", object$evalue[object$n]), "\n")
  cat("Max log(E):         ", sprintf("%.4f", max(object$log_evalue)), "\n")
  cat("Always-valid p:     ", sprintf("%.6f", 1 / max(exp(max(object$log_evalue)), 1)), "\n")
  if (object$rejected) {
    cat("Decision:           REJECT at observation ", object$rejection_time, "\n", sep = "")
  } else {
    cat("Decision:           Do not reject\n")
  }
  invisible(object)
}

#' Print an E-Value Monitor
#'
#' Display the current state of an e-value monitor, including the cumulative
#' sample size, current e-value, always-valid p-value, confidence sequence,
#' and rejection status.
#'
#' @param x An \code{"emonitor"} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.emonitor <- function(x, ...) {
  cat("E-value monitor\n")
  cat("  Current n:  ", x$n, "\n", sep = "")
  cat("  Current E:  ", sprintf("%.3f", exp(x$log_evalue)), "\n", sep = "")
  cat("  log(E):     ", sprintf("%.3f", x$log_evalue), "\n", sep = "")
  cat("  AV p-value: ", sprintf("%.6f", x$av_pvalue), "\n", sep = "")
  if (!is.null(x$cs)) {
    cat("  CS for delta: [", sprintf("%.4f", x$cs[1]),
        ", ", sprintf("%.4f", x$cs[2]), "]\n", sep = "")
  }
  if (x$rejected) {
    cat("  Status: REJECTED\n")
  } else {
    cat("  Status: Monitoring\n")
  }
  invisible(x)
}

#' Summarize an E-Value Monitor
#'
#' @param object An \code{"emonitor"} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{object}.
#'
#' @export
summary.emonitor <- function(object, ...) {
  print.emonitor(object, ...)
}

#' Print a Comparison Results Object
#'
#' Display a formatted table comparing the performance of different monitoring
#' methods, including Type I error rate, power, and average sample size under
#' both the null and alternative hypotheses.
#'
#' @param x An \code{"ecomparison"} object from \code{\link{simulate_comparison}}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.ecomparison <- function(x, ...) {
  cat("Monitoring method comparison\n")
  cat("============================\n")
  cat("Design: p_C =", x$design$p_C, ", p_T (alt) =", x$design$p_T_alt, "\n")
  cat("Nmax =", x$design$Nmax, "per arm,", x$design$n_looks, "looks, alpha =",
      x$design$alpha, "\n")
  cat("Replications:", x$design$nrep, "\n\n")

  res <- x$results
  fmt <- "  %-35s %8s %8s %8s %8s\n"
  cat(sprintf(fmt, "Method", "NullRej", "AltRej", "AvgN_0", "AvgN_1"))
  cat(sprintf(fmt, strrep("-", 35), strrep("-", 8), strrep("-", 8),
              strrep("-", 8), strrep("-", 8)))
  for (i in seq_len(nrow(res))) {
    cat(sprintf(fmt,
                res$method[i],
                sprintf("%.4f", res$null_rej[i]),
                sprintf("%.4f", res$alt_rej[i]),
                sprintf("%.1f", res$avg_n_null[i]),
                sprintf("%.1f", res$avg_n_alt[i])))
  }
  invisible(x)
}

#' Print a Confidence Sequence Object
#'
#' Display the current state of an always-valid confidence sequence,
#' including the number of observations and the latest confidence interval.
#'
#' @param x A \code{"confseq"} object from \code{\link{confseq_binary}}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.confseq <- function(x, ...) {
  cat("Confidence sequence for treatment effect\n")
  cat("  Observations: ", x$n, "\n", sep = "")
  cat("  Level:        ", 1 - x$alpha, "\n", sep = "")
  cat("  Current CS:   [", sprintf("%.4f", x$lower[x$n]),
      ", ", sprintf("%.4f", x$upper[x$n]), "]\n", sep = "")
  invisible(x)
}
