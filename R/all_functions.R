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
  cat("Components selected: nPC1 =", x$nPC_used$nPC1, ", nPC2 =", x$nPC_used$nPC2, "\n")
  
  if (!is.null(x$performance)) {
    cat("\nPerformance Metrics:\n")
    cat("  Exposure 1 - MISE:", round(x$performance$MISE1, 6), 
        ", Coverage:", round(x$performance$Coverage1, 3), "\n")
    cat("  Exposure 2 - MISE:", round(x$performance$MISE2, 6),
        ", Coverage:", round(x$performance$Coverage2, 3), "\n")
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
  if (!is.null(x$plots$p1) && !is.null(x$plots$p2)) {
    requireNamespace("gridExtra", quietly = TRUE)
    gridExtra::grid.arrange(x$plots$p1, x$plots$p2, ncol = 2)
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
  
  cat("\nExposure 1:\n")
  cat("  Components:", x$exposure1$nPC_used, "\n")
  if (!is.null(x$exposure1$performance)) {
    cat("  MSE:", round(x$exposure1$performance$MISE, 6), "\n")
    cat("  Coverage:", round(x$exposure1$performance$Coverage, 3), "\n")
  }
  
  if (!is.null(x$exposure2$performance$MISE)) {
    cat("\nExposure 2:\n")
    cat("  Components:", x$exposure2$nPC_used, "\n")
    cat("  MSE:", round(x$exposure2$performance$MISE, 6), "\n")
    cat("  Coverage:", round(x$exposure2$performance$Coverage, 3), "\n")
  }
  
  invisible(x)
}

#' @export
summary.mvfmr_separate <- function(object, ...) {
  print(object)
  
  cat("\nExposure 1 Coefficients:\n")
  print(round(object$exposure1$coefficients, 4))
  
  cat("\nExposure 2 Coefficients:\n")
  print(round(object$exposure2$coefficients, 4))
  
  invisible(object)
}

#' @export
plot.mvfmr_separate <- function(x, ...) {
  if (!is.null(x$plots$p1) && !is.null(x$plots$p2)) {
    requireNamespace("gridExtra", quietly = TRUE)
    gridExtra::grid.arrange(x$plots$p1, x$plots$p2, ncol = 2)
  } else if(!is.null(x$plots$p1)) {
    requireNamespace("gridExtra", quietly = TRUE)
    gridExtra::grid.arrange(x$plots$p1, ncol = 1)
  }
  invisible(NULL)
}

#' @export
coef.mvfmr_separate <- function(object, exposure = c(1, 2), ...) {
  exposure <- match.arg(as.character(exposure), c("1", "2"))
  
  if (exposure == "1") {
    return(object$exposure1$coefficients)
  } else {
    return(object$exposure2$coefficients)
  }
}

#' @export
vcov.mvfmr_separate <- function(object, exposure = c(1, 2), ...) {
  exposure <- match.arg(as.character(exposure), c("1", "2"))
  
  if (exposure == "1") {
    return(object$exposure1$vcov)
  } else {
    return(object$exposure2$vcov)
  }
}

