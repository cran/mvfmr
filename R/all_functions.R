# ============= FMVMR Package - S3 Methods and Helper Functions # =============

# ============= S3 METHODS FOR mvfmr # =============

#' @export
print.mvfmr <- function(x, ...) {
  cat("\nFunctional Multivariable MR Result\n")
  cat("==================================\n")
  cat("Exposures:", x$n_exposures, "\n")
  cat("Sample size:", x$n_observations, "\n")
  cat("Outcome:", x$outcome_type, "\n")
  cat("Method:", x$method, "\n")
  cat("Components selected:", paste(paste0("nPC", seq_len(x$n_exposures), " = ", x$nPC_used), collapse = ", "), "\n")

  if (!is.null(x$performance)) {
    cat("\nPerformance Metrics:\n")
    for (k in seq_len(x$n_exposures)) {
      cat("  Exposure", k, "- MISE:", round(x$performance$MISE[[k]], 6),
          ", Coverage:", round(x$performance$Coverage[[k]], 3), "\n")
    }
  }

  invisible(x)
}

#' @export
summary.mvfmr <- function(object, ...) {
  print(object)
  cat("\nCoefficients:\n")
  print(round(object$coefficients, 4))
  cat("\nUse plot() to visualize time-varying effects\n")
  invisible(object)
}

#' @export
plot.mvfmr <- function(x, ...) {
  if (!is.null(x$plots$effects) && length(x$plots$effects) > 0) {
    requireNamespace("gridExtra", quietly = TRUE)
    gridExtra::grid.arrange(grobs = x$plots$effects, ncol = ceiling(sqrt(length(x$plots$effects))))
  } else {
    message("Plots not available in this result object")
  }
  invisible(NULL)
}

#' @export
coef.mvfmr <- function(object, ...) {
  if (!is.null(object$coefficients)) {
    return(object$coefficients)
  }
  message("No coefficients available")
  return(NULL)
}

#' @export
vcov.mvfmr <- function(object, ...) {
  if (!is.null(object$vcov)) {
    return(object$vcov)
  }
  message("No variance-covariance matrix available")
  return(NULL)
}

# ============= S3 METHODS FOR fmvmr_separate # =============

#' @export
print.mvfmr_separate <- function(x, ...) {
  cat("\nSeparate Univariable MR Results\n")
  cat("================================\n")
  cat("Exposures:", x$n_exposures, "\n")
  cat("Separate instruments:", x$separate_instruments, "\n")
  cat("Outcome:", x$outcome_type, "\n")
  cat("Method:", x$method, "\n")

  for (k in seq_len(x$n_exposures)) {
    cat("\nExposure", k, ":\n")
    cat("  Components:", x$exposures[[k]]$nPC_used, "\n")
    if (!is.null(x$exposures[[k]]$performance)) {
      cat("  MSE:", round(x$exposures[[k]]$performance$MISE, 6), "\n")
      cat("  Coverage:", round(x$exposures[[k]]$performance$Coverage, 3), "\n")
    }
  }

  invisible(x)
}

#' @export
summary.mvfmr_separate <- function(object, ...) {
  print(object)

  for (k in seq_len(object$n_exposures)) {
    cat("\nExposure", k, "Coefficients:\n")
    print(round(object$exposures[[k]]$coefficients, 4))
  }

  invisible(object)
}

#' @export
plot.mvfmr_separate <- function(x, ...) {
  if (!is.null(x$plots$effects) && length(x$plots$effects) > 0) {
    requireNamespace("gridExtra", quietly = TRUE)
    gridExtra::grid.arrange(grobs = x$plots$effects, ncol = ceiling(sqrt(length(x$plots$effects))))
  }
  invisible(NULL)
}

#' @export
coef.mvfmr_separate <- function(object, exposure = 1, ...) {
  if (!exposure %in% seq_len(object$n_exposures)) {
    stop("`exposure` must be an integer between 1 and ", object$n_exposures)
  }

  object$exposures[[exposure]]$coefficients
}

#' @export
vcov.mvfmr_separate <- function(object, exposure = 1, ...) {
  if (!exposure %in% seq_len(object$n_exposures)) {
    stop("`exposure` must be an integer between 1 and ", object$n_exposures)
  }

  object$exposures[[exposure]]$vcov
}
