### Imports
library(MASS)
library(MplusAutomation)
library(lslx)
source("performancemeasures.R")


#' Valid Args
#'
#' Gives all combinations of arguments that constitute simulation settings used
#' in the study.
#'
#' @return A matrix with all possible combinations of arguments that can be
#' passed to 'get_simulation_arguments'.
#' @export

valid_args <- function() {
  combs <- data.frame(matrix(NA, ncol = 5, nrow = 40))
  colnames(combs) <- c("groups", "components", "differences", "CL", "PLdecrease")

  groups <- c(2, 4)
  components <- c(2, 4)
  differences <- c(4, 16)
  CL_PL <- list(c(0, 0), c(0.2, 0), c(0.4, 0), c(0, 0.2), c(0, 0.4))

  grd <- expand.grid(groups, components, differences, CL_PL)
  combs[, 1:3] <- grd[, 1:3]
  combs[, 4] <- unlist(grd[, 4])[c(TRUE, FALSE)]
  combs[, 5] <- unlist(grd[, 4])[c(FALSE, TRUE)]

  combs_sorted <- combs[with(combs, order(groups, components, differences, CL, PLdecrease)), ]
  row.names(combs_sorted) <- NULL

  return(combs_sorted)
}


#' Get Simulation Arguments
#'
#' Gets values for simulation arguments for a given simulation setting.
#'
#' @param groups Number of groups (2 or 4)
#' @param components Number of groups (2 or 4)
#' @param differences Number of loading differences (4 or 16)
#' @param CL Cross Loadings value (0, 0.2 or 0.4)
#' @param PLdecrease Primary Loading decrease value (0, 0.2 or 0.4)
#'
#' @return list with entries P1, P2, SIGMA1, SIGMA2; to be used for the simulation
#' @export

get_simulation_arguments <- function(groups, components, differences, CL, PLdecrease) {
  args <- list()

  if (groups == 2 & components == 2 & differences == 4 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, rep(sqrt(.6), 9), rep(0, 10))
    p2g1 <- c(sqrt(.6), rep(0, 9), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, rep(sqrt(.6), 7), rep(0, 10))
    p2g2 <- c(0, 0, sqrt(.6), rep(0, 7), rep(sqrt(.6), 10))
  } else if (groups == 2 & components == 2 & differences == 4 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), .2, rep(0, 9))
    p2g1 <- c(.2, rep(0, 9), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), .2, rep(0, 5))
    p2g2 <- c(rep(0, 4), .2, rep(0, 5), rep(sqrt(.6), 10))
  } else if (groups == 2 & components == 2 & differences == 4 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), .4, rep(0, 9))
    p2g1 <- c(.4, rep(0, 9), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), .4, rep(0, 5))
    p2g2 <- c(rep(0, 4), .4, rep(0, 5), rep(sqrt(.6), 10))
  } else if (groups == 2 & components == 2 & differences == 16 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, 0, rep(sqrt(.6), 10), rep(0, 8))
    p2g1 <- c(sqrt(.6), sqrt(.6), rep(0, 10), rep(sqrt(.6), 8))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, 0, rep(sqrt(.6), 6), 0, 0, rep(sqrt(.6), 2), rep(0, 6))
    p2g2 <- c(0, 0, rep(sqrt(.6), 2), rep(0, 6), rep(sqrt(.6), 2), 0, 0, rep(sqrt(.6), 6))
  } else if (groups == 2 & components == 2 & differences == 16 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), rep(.2, 4), rep(0, 6))
    p2g1 <- c(rep(.2, 4), rep(0, 6), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), rep(.2, 4), 0, 0)
    p2g2 <- c(rep(0, 4), rep(.2, 4), 0, 0, rep(sqrt(.6), 10))
  } else if (groups == 2 & components == 2 & differences == 16 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), rep(.4, 4), rep(0, 6))
    p2g1 <- c(rep(.4, 4), rep(0, 6), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), rep(.4, 4), 0, 0)
    p2g2 <- c(rep(0, 4), rep(.4, 4), 0, 0, rep(sqrt(.6), 10))
  } else if (groups == 2 & components == 2 & differences == 4 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(sqrt(.6) - .2, rep(sqrt(.6), 9), rep(0, 10))
    p2g1 <- c(rep(0, 10), sqrt(.6) - .2, rep(sqrt(.6), 9))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), sqrt(.6) - .2, rep(sqrt(.6), 5), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), sqrt(.6) - .2, rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 2 & differences == 4 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(sqrt(.6) - .4, rep(sqrt(.6), 9), rep(0, 10))
    p2g1 <- c(rep(0, 10), sqrt(.6) - .4, rep(sqrt(.6), 9))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), sqrt(.6) - .4, rep(sqrt(.6), 5), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), sqrt(.6) - .4, rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 2 & differences == 16 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 6), rep(0, 10))
    p2g1 <- c(rep(0, 10), rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 6))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 2), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 2))
  } else if (groups == 2 & components == 2 & differences == 16 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 6), rep(0, 10))
    p2g1 <- c(rep(0, 10), rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 6))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 2), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 2))
  } else if (groups == 2 & components == 4 & differences == 4 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, rep(sqrt(.6), 4), rep(0, 15))
    p2g1 <- c(sqrt(.6), rep(0, 4), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, rep(sqrt(.6), 2), rep(0, 15))
    p2g2 <- c(rep(0, 2), sqrt(.6), 0, 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 4 & differences == 4 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), .2, rep(0, 14))
    p2g1 <- c(.2, rep(0, 4), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), 0, 0, .2, rep(0, 12))
    p2g2 <- c(rep(0, 2), .2, 0, 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 4 & differences == 4 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), .4, rep(0, 14))
    p2g1 <- c(.4, rep(0, 4), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), 0, 0, .4, rep(0, 12))
    p2g2 <- c(rep(0, 2), .4, 0, 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 4 & differences == 16 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, rep(sqrt(.6), 5), rep(0, 14))
    p2g1 <- c(sqrt(.6), rep(0, 5), rep(sqrt(.6), 4), rep(0, 10))
    p3g1 <- c(rep(0, 11), rep(sqrt(.6), 5), rep(0, 4))
    p4g1 <- c(rep(0, 10), sqrt(.6), rep(0, 5), rep(sqrt(.6), 4))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, rep(sqrt(.6), 2), 0, sqrt(.6), rep(0, 13))
    p2g2 <- c(rep(0, 2), sqrt(.6), 0, 0, sqrt(.6), 0, rep(sqrt(.6), 3), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 2), 0, rep(sqrt(.6), 2), 0, sqrt(.6), rep(0, 3))
    p4g2 <- c(rep(0, 12), sqrt(.6), 0, 0, sqrt(.6), 0, rep(sqrt(.6), 3))
  } else if (groups == 2 & components == 4 & differences == 16 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), rep(.2, 2), rep(0, 13))
    p2g1 <- c(rep(.2, 2), rep(0, 3), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(.2, 2), rep(0, 3))
    p4g1 <- c(rep(0, 10), rep(.2, 2), rep(0, 3), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), 0, 0, rep(.2, 2), rep(0, 11))
    p2g2 <- c(rep(0, 2), rep(.2, 2), 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 2), rep(.2, 2), 0)
    p4g2 <- c(rep(0, 12), rep(.2, 2), 0, rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 4 & differences == 16 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), rep(.4, 2), rep(0, 13))
    p2g1 <- c(rep(.4, 2), rep(0, 3), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(.4, 2), rep(0, 3))
    p4g1 <- c(rep(0, 10), rep(.4, 2), rep(0, 3), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), 0, 0, rep(.4, 2), rep(0, 11))
    p2g2 <- c(rep(0, 2), rep(.4, 2), 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 2), rep(.4, 2), 0)
    p4g2 <- c(rep(0, 12), rep(.4, 2), 0, rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 4 & differences == 4 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(sqrt(.6) - .2, rep(sqrt(.6), 4), rep(0, 15))
    p2g1 <- c(rep(0, 5), sqrt(.6) - .2, rep(sqrt(.6), 4), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), sqrt(.6) - .2, rep(sqrt(.6), 2), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), sqrt(.6) - .2, rep(sqrt(.6), 2), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 4 & differences == 4 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(sqrt(.6) - .4, rep(sqrt(.6), 4), rep(0, 15))
    p2g1 <- c(rep(0, 5), sqrt(.6) - .4, rep(sqrt(.6), 4), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), sqrt(.6) - .4, rep(sqrt(.6), 2), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), sqrt(.6) - .4, rep(sqrt(.6), 2), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 2 & components == 4 & differences == 16 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3), rep(0, 15))
    p2g1 <- c(rep(0, 5), rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6))
  } else if (groups == 2 & components == 4 & differences == 16 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3), rep(0, 15))
    p2g1 <- c(rep(0, 5), rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6))
  } else if (groups == 4 & components == 2 & differences == 4 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, rep(sqrt(.6), 9), rep(0, 10))
    p2g1 <- c(sqrt(.6), rep(0, 9), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, rep(sqrt(.6), 7), rep(0, 10))
    p2g2 <- c(0, 0, sqrt(.6), rep(0, 7), rep(sqrt(.6), 10))
  } else if (groups == 4 & components == 2 & differences == 4 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), .2, rep(0, 9))
    p2g1 <- c(.2, rep(0, 9), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), .2, rep(0, 5))
    p2g2 <- c(rep(0, 4), .2, rep(0, 5), rep(sqrt(.6), 10))
  } else if (groups == 4 & components == 2 & differences == 4 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), .4, rep(0, 9))
    p2g1 <- c(.4, rep(0, 9), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), .4, rep(0, 5))
    p2g2 <- c(rep(0, 4), .4, rep(0, 5), rep(sqrt(.6), 10))
  } else if (groups == 4 & components == 2 & differences == 16 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, 0, rep(sqrt(.6), 10), rep(0, 8))
    p2g1 <- c(sqrt(.6), sqrt(.6), rep(0, 10), rep(sqrt(.6), 8))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, 0, rep(sqrt(.6), 6), 0, 0, rep(sqrt(.6), 2), rep(0, 6))
    p2g2 <- c(0, 0, rep(sqrt(.6), 2), rep(0, 6), rep(sqrt(.6), 2), 0, 0, rep(sqrt(.6), 6))
  } else if (groups == 4 & components == 2 & differences == 16 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), rep(.2, 4), rep(0, 6))
    p2g1 <- c(rep(.2, 4), rep(0, 6), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), rep(.2, 4), 0, 0)
    p2g2 <- c(rep(0, 4), rep(.2, 4), 0, 0, rep(sqrt(.6), 10))
  } else if (groups == 4 & components == 2 & differences == 16 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 10), rep(.4, 4), rep(0, 6))
    p2g1 <- c(rep(.4, 4), rep(0, 6), rep(sqrt(.6), 10))

    # P2
    p1g2 <- c(rep(sqrt(.6), 10), rep(0, 4), rep(.4, 4), 0, 0)
    p2g2 <- c(rep(0, 4), rep(.4, 4), 0, 0, rep(sqrt(.6), 10))
  } else if (groups == 4 & components == 2 & differences == 4 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(sqrt(.6) - .2, rep(sqrt(.6), 9), rep(0, 10))
    p2g1 <- c(rep(0, 10), sqrt(.6) - .2, rep(sqrt(.6), 9))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), sqrt(.6) - .2, rep(sqrt(.6), 5), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), sqrt(.6) - .2, rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 2 & differences == 4 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(sqrt(.6) - .4, rep(sqrt(.6), 9), rep(0, 10))
    p2g1 <- c(rep(0, 10), sqrt(.6) - .4, rep(sqrt(.6), 9))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), sqrt(.6) - .4, rep(sqrt(.6), 5), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), sqrt(.6) - .4, rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 2 & differences == 16 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 6), rep(0, 10))
    p2g1 <- c(rep(0, 10), rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 6))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 2), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), rep(sqrt(.6) - .2, 4), rep(sqrt(.6), 2))
  } else if (groups == 4 & components == 2 & differences == 16 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 6), rep(0, 10))
    p2g1 <- c(rep(0, 10), rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 6))

    # P2
    p1g2 <- c(rep(sqrt(.6), 4), rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 2), rep(0, 10))
    p2g2 <- c(rep(0, 10), rep(sqrt(.6), 4), rep(sqrt(.6) - .4, 4), rep(sqrt(.6), 2))
  } else if (groups == 4 & components == 4 & differences == 4 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, rep(sqrt(.6), 4), rep(0, 15))
    p2g1 <- c(sqrt(.6), rep(0, 4), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, rep(sqrt(.6), 2), rep(0, 15))
    p2g2 <- c(rep(0, 2), sqrt(.6), 0, 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 4 & differences == 4 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), .2, rep(0, 14))
    p2g1 <- c(.2, rep(0, 4), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), 0, 0, .2, rep(0, 12))
    p2g2 <- c(rep(0, 2), .2, 0, 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 4 & differences == 4 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), .4, rep(0, 14))
    p2g1 <- c(.4, rep(0, 4), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), 0, 0, .4, rep(0, 12))
    p2g2 <- c(rep(0, 2), .4, 0, 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 4 & differences == 16 & CL == 0 & PLdecrease == 0) {
    # P1
    p1g1 <- c(0, rep(sqrt(.6), 5), rep(0, 14))
    p2g1 <- c(sqrt(.6), rep(0, 5), rep(sqrt(.6), 4), rep(0, 10))
    p3g1 <- c(rep(0, 11), rep(sqrt(.6), 5), rep(0, 4))
    p4g1 <- c(rep(0, 10), sqrt(.6), rep(0, 5), rep(sqrt(.6), 4))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), 0, rep(sqrt(.6), 2), 0, sqrt(.6), rep(0, 13))
    p2g2 <- c(rep(0, 2), sqrt(.6), 0, 0, sqrt(.6), 0, rep(sqrt(.6), 3), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 2), 0, rep(sqrt(.6), 2), 0, sqrt(.6), rep(0, 3))
    p4g2 <- c(rep(0, 12), sqrt(.6), 0, 0, sqrt(.6), 0, rep(sqrt(.6), 3))
  } else if (groups == 4 & components == 4 & differences == 16 & CL == 0.2 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), rep(.2, 2), rep(0, 13))
    p2g1 <- c(rep(.2, 2), rep(0, 3), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(.2, 2), rep(0, 3))
    p4g1 <- c(rep(0, 10), rep(.2, 2), rep(0, 3), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), rep(0, 2), rep(.2, 2), rep(0, 11))
    p2g2 <- c(rep(0, 2), rep(.2, 2), 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 2), rep(.2, 2), 0)
    p4g2 <- c(rep(0, 12), rep(.2, 2), 0, rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 4 & differences == 16 & CL == 0.4 & PLdecrease == 0) {
    # P1
    p1g1 <- c(rep(sqrt(.6), 5), rep(.4, 2), rep(0, 13))
    p2g1 <- c(rep(.4, 2), rep(0, 3), rep(sqrt(.6), 5), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(.4, 2), rep(0, 3))
    p4g1 <- c(rep(0, 10), rep(.4, 2), rep(0, 3), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 5), rep(0, 2), rep(.4, 2), rep(0, 11))
    p2g2 <- c(rep(0, 2), rep(.4, 2), 0, rep(sqrt(.6), 5), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 2), rep(.4, 2), 0)
    p4g2 <- c(rep(0, 12), rep(.4, 2), 0, rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 4 & differences == 4 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(sqrt(.6) - .2, rep(sqrt(.6), 4), rep(0, 15))
    p2g1 <- c(rep(0, 5), sqrt(.6) - .2, rep(sqrt(.6), 4), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), sqrt(.6) - .2, rep(sqrt(.6), 2), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), sqrt(.6) - .2, rep(sqrt(.6), 2), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 4 & differences == 4 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(sqrt(.6) - .4, rep(sqrt(.6), 4), rep(0, 15))
    p2g1 <- c(rep(0, 5), sqrt(.6) - .4, rep(sqrt(.6), 4), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6), 5))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), sqrt(.6) - .4, rep(sqrt(.6), 2), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), sqrt(.6) - .4, rep(sqrt(.6), 2), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 5), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 5))
  } else if (groups == 4 & components == 4 & differences == 16 & CL == 0 & PLdecrease == 0.2) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3), rep(0, 15))
    p2g1 <- c(rep(0, 5), rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6) - .2, 2), rep(sqrt(.6), 3))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 2), rep(sqrt(.6) - .2, 2), sqrt(.6))
  } else if (groups == 4 & components == 4 & differences == 16 & CL == 0 & PLdecrease == 0.4) {
    # P1
    p1g1 <- c(rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3), rep(0, 15))
    p2g1 <- c(rep(0, 5), rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3), rep(0, 10))
    p3g1 <- c(rep(0, 10), rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3), rep(0, 5))
    p4g1 <- c(rep(0, 15), rep(sqrt(.6) - .4, 2), rep(sqrt(.6), 3))

    # P2
    p1g2 <- c(rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6), rep(0, 15))
    p2g2 <- c(rep(0, 5), rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6), rep(0, 10))
    p3g2 <- c(rep(0, 10), rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6), rep(0, 5))
    p4g2 <- c(rep(0, 15), rep(sqrt(.6), 2), rep(sqrt(.6) - .4, 2), sqrt(.6))
  } else {
    stop("Invalid combination of arguments. Call 'valid_args()' to see all
         valid combinations of arguments for simulation.")
  }


  if (components == 2) {
    P1 <- cbind(p1g1, p2g1)
    P2 <- cbind(p1g2, p2g2)
  } else if (components == 4) {
    P1 <- cbind(p1g1, p2g1, p3g1, p4g1)
    P2 <- cbind(p1g2, p2g2, p3g2, p4g2)
  }


  PSI1 <- diag(1 - rowSums(P1^2))
  SIGMA1 <- P1 %*% t(P1) + PSI1
  PSI2 <- diag(1 - rowSums(P2^2))
  SIGMA2 <- P2 %*% t(P2) + PSI2

  args$P1 <- P1
  args$P2 <- P2
  args$SIGMA1 <- SIGMA1
  args$SIGMA2 <- SIGMA2

  return(args)
}

#' Simulation
#'
#' Runs one simulation for a given setting.
#'
#' @param groups Number of groups (2 or 4)
#' @param components Number of groups (2 or 4)
#' @param differences Number of loading differences (4 or 16)
#' @param CL Cross Loadings value (0, 0.2 or 0.4)
#' @param PLdecrease Primary Loading decrease value (0, 0.2 or 0.4)
#' @param Illustrative Logical: single illustrative example or full simulation
#' @param seed Seed to use (default 45468)
#'
#' @return Results for the given simulation setting.
#' @export

simulation_LSLX <- function(groups, components, differences, CL, PLdecrease, Illustrative = F, seed = 45468) {
  N <- c(50, 200, 600, 1000)
  if (Illustrative){
    final_result <- list()  
  } else {
  final_result <- data.frame()
  }
  set.seed(seed)

  args <- get_simulation_arguments(groups, components, differences, CL, PLdecrease)
  P1 <- args$P1
  P2 <- args$P2
  SIGMA1 <- args$SIGMA1
  SIGMA2 <- args$SIGMA2


  for (i in 1:length(N)) {
    getValue <- function() {
      
      #make data in lslx format
      if (groups == 2) {
        data1 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data2 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
        data_lslx <- data.frame(rbind(data1,data2))
        grouplabel <- as.factor(rep(c("G1","G2"),c(N[i],N[i])))
        data_lslx <- data.frame(cbind(data_lslx,grouplabel))
      } else if (groups == 4) {
        data1 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data2 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
        data3 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data4 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
        data_lslx <- data.frame(rbind(data1,data2,data3,data4))
        grouplabel <- as.factor(rep(c("G1","G2","G3","G4"),c(rep(N[i],groups))))
        data_lslx <- data.frame(cbind(data_lslx,grouplabel))
      }
      
      #define lslx model
      if (components==2){
        samplex <- seq(1:20)
        #line1 <-paste("f1:=> ",paste("X",samplex[1:10],collapse = "+",sep = ""))
        #line2 <-paste("f2:=> ",paste("X",samplex[11:20],collapse = "+",sep = ""))
        line3 <-paste("f1:~> ",paste("X",samplex[1:20],collapse = "+",sep = ""))
        line4 <-paste("f2:~> ",paste("X",samplex[1:20],collapse = "+",sep = ""))
        linex <- paste("f1  <=>  f2")
        line5 <-paste("f1  <=> fix(1)* f1")
        line6 <-paste("f2  <=> fix(1)* f2")
        model_mgfa <- paste(paste(line3,line4,linex,line5,line6,collapse = "\n", sep="\n"),collapse = "\n","\n")
      } else {
        samplex <- seq(1:20)
        #line1 <-paste("f1:=> ",paste("X",samplex[1:5],collapse = "+",sep = ""))
        #line2 <-paste("f2:=> ",paste("X",samplex[6:10],collapse = "+",sep = ""))
        #line3 <-paste("f3:=> ",paste("X",samplex[11:15],collapse = "+",sep = ""))
        #line4 <-paste("f4:=> ",paste("X",samplex[16:20],collapse = "+",sep = ""))       
        line5 <-paste("f1:~> ",paste("X",samplex[1:20],collapse = "+",sep = ""))
        line6 <-paste("f2:~> ",paste("X",samplex[c(1:10,11:20)],collapse = "+",sep = ""))
        line7 <-paste("f3:~> ",paste("X",samplex[c(1:15,16:20)],collapse = "+",sep = ""))
        line8 <-paste("f4:~> ",paste("X",samplex[1:20],collapse = "+",sep = ""))
        linex1 <- paste("f1  <=>  f2")
        linex2 <- paste("f1  <=>  f3")
        linex3 <- paste("f3  <=>  f2")
        line9 <-paste("f1  <=> fix(1)* f1")
        line10 <-paste("f2  <=> fix(1)* f2")
        line11 <-paste("f3  <=> fix(1)* f3")
        line12 <-paste("f4  <=> fix(1)* f4")
        model_mgfa <- paste(paste(line5,line6,line7,line8,linex1,linex2,linex3,line9,line10,line11,line12,collapse = "\n", sep="\n"),collapse = "\n","\n")
      }
      
      #run lxlx analysis
      lslx_mgfa <- lslx$new(model = model_mgfa,
                            data = data_lslx,
                            group_variable = "grouplabel",
                            reference_group = "G1")#,
      if (groups==4){
        lslx_mgfa$penalize_heterogeneity(
          block = c("y<-f","y<-1"),
          group = "G2"
        )
        lslx_mgfa$penalize_heterogeneity(
          block = c("y<-f","y<-1"),
          group = "G3"
        )
        lslx_mgfa$penalize_heterogeneity(
          block = c("y<-f","y<-1"),
          group = "G4"
        )
      } else {
        lslx_mgfa$penalize_heterogeneity(
          block = c("y<-f","y<-1"),
          group = "G2"
        )
      }
      #lslx_mgfa$free_block(block = "f<-1")
      lslx_mgfa$fit(penalty_method = "mcp", 
                    lambda_grid = seq(.01, .60, .01), 
                    delta_grid = c(1.5, 3.0, Inf))
      #lslx_mgfa$summarize(selector = "bic", interval = FALSE, include_faulty = T)

        error_count<-0
        result <- tryCatch(
          {
            lslx_mgfa$summarize(selector = "bic")   # or whichever call triggers the error
          },
          error = function(e) {
            if (grepl("PL estimate under EACH penalty level.*nonconverged result",
                      conditionMessage(e),
                      ignore.case = TRUE)) {
              error_count <<- 1
              #message("Detected the nonconverged-all-lambda error.")
              return(NULL)
            }
            if (grepl("Approximated Hessian is not convex under EVERY convexity level",
                      conditionMessage(e),
                      ignore.case = TRUE)) {
              error_count <<- 1
              #message("Detected the nonconverged-all-lambda error.")
              return(NULL)
            }
          }
        )
        if (error_count==0){
          lslxcoefs <- lslx_mgfa$extract_coefficient(selector = "bic", include_faulty = T)
          #POSTPROCESS lslx OUTPUT
          loadings <- list()
          for (g in 1:groups){
            loadings[[g]] <- P1
            for (r in 1:components){
              strRG <- paste("<-f",r,"/",levels(grouplabel)[g],sep="")
              index <- grep(strRG,names(lslxcoefs))
              if (g==1){#reference group
                loadings[[g]][,r] <- lslxcoefs[index]
              } else {
                loadings[[g]][,r] <- lslxcoefs[index]+loadings[[1]][,r] 
              }
            }
          }
          lslxmodel <- list()
          lslxmodel$loadings <- loadings
          if (Illustrative){
            return(lslxmodel)
          } else {
            #calculation performance measures
            parameters <- performancemeasures(lslxmodel,P1,P2,components,groups)[[1]]    
          }
        } else {
        parameters <- c(rep(1,6),0)    
        }
          return(c(parameters,error_count))
    }
    if (Illustrative){
      result <- getValue()
      tmp <- list('loadings',result, 'n', N[i], 'nrpcs', components, 
                           'nrgroups', groups)
      final_result[[i]] <- tmp
    } else {
      result <- as.data.frame(t(replicate(10, getValue())))
      result$size <- c(rep(N[i], 10))
      result$nrpcs <- components
      result$nrgroups <- groups
      final_result <- rbind(final_result, result)
      save(final_result, file = paste("LSLX_result25.RData",sep=''))
    }
  }
  return(final_result)
}

#' Run Simulations
#'
#' Runs simulations for all simulation settings.
#' @param Illustrative Logical to indicate if single illustrative example or full 
#' simulation is wanted
#'
#' @return A list of results, one per simulation setting.
#' @export

run_simulations <- function() {
  results <- list()
  combs <- valid_args()
  combs_nr <- nrow(combs)

  for (i in 25:25) {
    comb <- combs[i, ]
    results[[i]] <- do.call(simulation_LSLX, as.list(as.numeric(comb)))
    ## save intermediary results
    res_lslx <- results[[i]]
    save(res_lslx, file = paste("LSLX_result",i,".RData",sep=''))
    ##
  }
  ##
  ##
  return(results)
}
LSLXresuls <- run_simulations()
save(LSLXresults,"LSLXresultsLoadings.RData")

#' Join Simulation Results
#' 
#' Organises Simulation Results.
#'
#' @param results A list of results, generated by 'run_simulations'
#'
#' @return A dataframe with all results joined.
#' @export

join_simulation_results <- function(results) {
  result <- c()
  G2R2_ids <- c(1:10) # 2 groups 2 components positions in results
  G2R4_ids <- c(11:20) # 2 groups 4 components positions in results
  G4R2_ids <- c(21:30) # 4 groups 2 components positions in results
  G4R4_ids <- c(31:40) # 4 groups 4 components positions in results
  
  # one big dataframe for 2 groups 2 components
  dataG2R2 <- c()
  for (i in G2R2_ids){
    load(paste("LSLX_result",i,".RData",sep=''))
    dataG2R2 <- rbind(dataG2R2,res_lslx)
  }
  # one big dataframe for 2 groups 4 components
  dataG2R4 <- c()
  for (i in G2R4_ids){
    load(paste("LSLX_result",i,".RData",sep=''))
    dataG2R4 <- rbind(dataG2R4,res_lslx)
  }
  # one big dataframe for 4 groups 2 components
  dataG4R2 <- c()
  for (i in G4R2_ids){
    load(paste("LSLX_result",i,".RData",sep=''))
    dataG4R2 <- rbind(dataG4R2,res_lslx)
  }
  # one big dataframe for 4 groups 4 components
  dataG4R4 <- c()
  for (i in G4R4_ids){
    load(paste("LSLX_result",i,".RData",sep=''))
    dataG4R4 <- rbind(dataG4R4,res_lslx)
  }
  
  # join
  result <- rbind(dataG2R2, dataG2R4, dataG4R2, dataG4R4)
  #get loading pattern conditions
  combs <- valid_args()
  nr_n <- 4
  nrreplics <- rep(c(50*nr_n,10*nr_n),c(4,36))
  CLpattern <- rep(combs$CL,nrreplics)
  PLpattern <- rep(combs$PLdecrease,nrreplics)
  nrdiff <- rep(combs$differences,nrreplics)
  result <- cbind(result,CLpattern,PLpattern,nrdiff)
  
  # rename columns
  labels <- c(
    "FP_PC1", "FP_PC2", "FN_PC1", "FN_PC2",
    "FUSION_FP", "FUSION_FN", "Tucker","Converg",
    "n", "R", "G", "CL", "PL","NrDiff"
  )
  colnames(result) <- labels
  
  # Format output
  output <- list()
  output$G2R2 <- dataG2R2
  output$G2R4 <- dataG2R4
  output$G4R2 <- dataG4R2
  output$G4R4 <- dataG4R4
  output$RESULT <- result
  
  return(output)
}


###scripts for producing summary table
#make table to compare with De Roover's results
RESULT <- join_simulation_results(c())
RESULT <- RESULT$RESULT
TABLE <- matrix(nrow=16,ncol=3)  #%datasets no FP/FN for fusion and simple structure, everything correct
colnames(TABLE) <- c('FUSION', 'SIMPLE S.', 'BOTH')

TABLE_percondition <- matrix(0,nrow=40,ncol=10)
#"n", "R", "G", "CL", "PL","NrDiff"
colnames(TABLE_percondition) <- c('FUSION', 'SIMPLE S.', 'BOTH',"TUCKER","N","R","G","CL","PL","NrDiff")
for (n in c(50,200,600,1000)){
  RESULT_n <- RESULT[RESULT$n==n,]
  U <- 0
  nrreplics <- 50
  for (c in 1:40){
  L <- U+1
  U <- L+nrreplics-1
  TABLE_percondition[c,5:10] <- as.matrix(RESULT_n[L,8:13])
  TABLE_percondition[c,1] <- sum(rowSums(RESULT_n[L:U,5:6]==0)==2)#ISF
  TABLE_percondition[c,2] <- sum(rowSums(RESULT_n[L:U,1:4]==0)==4)#ISF
  TABLE_percondition[c,3] <- sum(rowSums(RESULT_n[L:U,1:6]==0)==6)#ISF
}
}
TABLE_percondition[TABLE_percondition[,9]==0.4,]

#find number of datasets 100%correct: OVERALL
fusion_ISF <- sum(rowSums(RESULT[,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[,1:6]==0)==6)#ISF
TABLE[16,] <- c(fusion_ISF,simple_ISF,both_ISF)/2240

#find number of datasets 100%correct: Per nr of groups
Rvec <- c(2,4)
Gvec <- c(2,4)
for (g in 1:2){
  sel <- which(RESULT$G==Gvec[g])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[g,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
}
#find number of datasets 100%correct: Per sample size
Nvec <- c(50,200,600,1000)
for (n in 1:4){
  sel <- which(RESULT$n==Nvec[n])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[2+n,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
}
#find number of datasets 100%correct: Per group size
for (q in 1:2){
  sel <- which(RESULT$R==Rvec[q])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[6+q,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
}
#specific conditions
#1.Primary loading shifts: 4800 data sets
sel <- which(RESULT$CL==0 | RESULT$PL==0)
fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
TABLE[9,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
#2.remaining conditions
conds <- c('CL0.4','CL0.2','PL0.4','PL0.2')
sel <- which(RESULT$CL==0.4)
fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
TABLE[10,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
sel <- which(RESULT$CL==0.2)
fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
TABLE[11,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
sel <- which(RESULT$PL==0.4)
fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
TABLE[12,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
sel <- which(RESULT$PL==0.2)
fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
TABLE[13,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
#find number of datasets 100%correct: Per nr of between group differences
Dvec <- c(4,16)
for (d in 1:2){
  sel <- which(RESULT$NrDiff==Dvec[d])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[13+d,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
}
write.table(TABLE, file='LSLXsimresults.txt',sep = '\t',dec = ",")

