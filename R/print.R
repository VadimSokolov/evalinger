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

#' @export
summary.emonitor <- function(object, ...) {
  print.emonitor(object, ...)
}

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

#' @export
print.confseq <- function(x, ...) {
  cat("Confidence sequence for treatment effect\n")
  cat("  Observations: ", x$n, "\n", sep = "")
  cat("  Level:        ", 1 - x$alpha, "\n", sep = "")
  cat("  Current CS:   [", sprintf("%.4f", x$lower[x$n]),
      ", ", sprintf("%.4f", x$upper[x$n]), "]\n", sep = "")
  invisible(x)
}
