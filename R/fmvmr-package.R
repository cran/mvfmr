#' mvfmr: Multivariable Functional Mendelian Randomization
#'
#' @description
#' Implements Multivariable Functional Mendelian randomization to estimate
#' time-varying causal effects of multiple correlated longitudinal exposures.
#'
#' @name mvfmr-package
#' @aliases mvfmr-package _PACKAGE
#'
#' @author Nicole Fontana
#'
#' @import fdapace
#' @import ggplot2
#' @importFrom stats coef vcov lm pnorm qnorm pchisq quantile
#' @importFrom stats rbinom rnorm runif residuals sd var
#' @importFrom stats na.omit fitted nlminb predict spline time
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom progress progress_bar
#' @importFrom pROC roc auc
#' @importFrom glmnet cv.glmnet
#' @importFrom gridExtra grid.arrange
NULL
