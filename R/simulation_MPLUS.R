### Imports
library(MASS)
library(MplusAutomation)
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
#' @param seed Seed to use (default 45468)
#'
#' @return Results for the given simulation setting.
#' @export

simulation_MPLUS <- function(groups, components, differences, CL, PLdecrease, Illustrative=T, seed = 45468) {
  N <- c(50, 200, 600, 1000)
  final_result <- data.frame()
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
      if (groups == 2) {
        data1 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data2 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
      } else if (groups == 4) {
        data1 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data2 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
        data3 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data4 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
      }

      G <- groups
      R <- components

      DATA <- cbind(data.frame(rbind(data1,data2)),rep(c(1:G),rep(N[i],G)))
      colnames(DATA) <- c("v1","v2","v3","v4","v5","v6","v7","v8","v9","v10",
                          "v11","v12","v13","v14","v15","v16","v17","v18","v19","v20","G")
      
      
      #MPLUS
      if (G==2){
        varstr <- paste("usevariables = v1 - v20; GROUPING = G(1=group1 2=group2);")
      } else {
        varstr <- paste("usevariables = v1 - v20; GROUPING = G(1=group1 2=group2 3=group3 4=group4);")
      }
      if (R==2){
        modelstr <- paste("f1-f2 BY v1-v20 (*1);")
      } else {
        modelstr <- paste("f1-f4 BY v1-v20 (*1);")
      }
      m.penefa <- mplusObject(
        TITLE = "Exploratory Penalized Factor Analysis",
        VARIABLE = varstr,#"usevariables = v1 - v20; GROUPING = G(1=group1 2=group2);",
        ANALYSIS = "alignment = free;",
        MODEL = modelstr,
        OUTPUT = "align;",
        #      usevariables = colnames(d),
        rdata = DATA)
      
      m.penefa.fit <- mplusModeler(m.penefa, dataout="2gr4pcs.dat", 
                                   modelout="2gr4pcs.inp",run = TRUE, 
                                   writeData = "ifmissing")
      
      #POSTPROCESS MPLUS OUTPUT
      loadingsparams <- m.penefa.fit$results$parameters$unstandardized
      loadings <- list()
      for (g in 1:G){
        loadings[[g]] <- P1
        for (r in 1:R){
          strR <- paste("F",r,".BY",sep="")
          strG <- paste("GROUP",g,sep="")
          if (N[i] < prod(dim(P1))){
            index <- loadingsparams$paramHeader==strR & loadingsparams$Group==strG
            loadings[[g]][,r] <- as.numeric(loadingsparams$est[index])
          } else {
            index <- loadingsparams$paramHeader==strR & loadingsparams$Group==strG
            loadings[[g]][,r] <- as.numeric(loadingsparams$est[index])
            if (dim(loadingsparams)[2]>4){
              loadings[[g]][loadingsparams$pval[index] > 0.05,r] <- 0
            }
          }
        }
      }
      mplusmodel <- list()
      mplusmodel$loadings <- loadings
      
      #calculation performance measures
      parameters <- performancemeasures(mplusmodel,P1,P2,R,G)[[1]]    
      #round(cbind(P1,selmodel$loadings[[1]],P2,selmodel$loadings[[2]]),2)
      if (Illustrative){
        return(loadings)
      } else {
      return(parameters)
      }
    }
    if (Illustrative){
      result <- getValue()
      tmp <- list('loadings',result, 'n', N[i], 'nrpcs', components, 
                  'nrgroups', groups)
      final_result[[i]] <- tmp
    } else {
      result <- as.data.frame(t(replicate(50, getValue())))
      result$size <- c(rep(N[i], 50))
      result$nrpcs <- components
      result$nrgroups <- groups
      final_result <- rbind(final_result, result)
    }
  }
  return(final_result)
}

#' Run Simulations
#'
#' Runs simulations for all simulation settings.
#'
#' @return A list of results, one per simulation setting.
#' @export

run_simulations <- function() {
  results <- list()
  combs <- valid_args()
  combs_nr <- nrow(combs)

  for (i in 1:combs_nr) {
    comb <- combs[i, ]
    results[[i]] <- do.call(simulation_MPLUS, as.list(as.numeric(comb)))
    ## save intermediary results
    res_mplus <- results[[i]]
    save(res_mplus, file = paste("MPLUS_resultLoadings",i,".RData",sep=''))
    ##
  }
  ##
  ##
  return(results)
}
MPLUSresuls <- run_simulations()
save(results,"MPLUSresultsLoadings.RData")

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
    dataG2R2 <- rbind(dataG2R2,results[[i]]$res_mplus)
  }
  # one big dataframe for 2 groups 4 components
  dataG2R4 <- c()
  for (i in G2R4_ids){
    dataG2R4 <- rbind(dataG2R4,results[[i]]$res_mplus)
  }
  # one big dataframe for 4 groups 2 components
  dataG4R2 <- c()
  for (i in G4R2_ids){
    dataG4R2 <- rbind(dataG4R2,results[[i]]$res_mplus)
  }
  # one big dataframe for 4 groups 4 components
  dataG4R4 <- c()
  for (i in G4R4_ids){
    dataG4R4 <- rbind(dataG4R4,results[[i]]$res_mplus)
  }

  # join
  result <- rbind(dataG2R2, dataG2R4, dataG4R2, dataG4R4)
  #get loading pattern conditions
  combs <- valid_args()
  nrreplics <- 50
  nr_n <- 4
  CLpattern <- rep(combs$CL,c(rep(nrreplics*nr_n,40)))
  PLpattern <- rep(combs$PLdecrease,c(rep(nrreplics*nr_n,40)))
  nrdiff <- rep(combs$differences,c(rep(nrreplics*nr_n,40)))
  result <- cbind(result,CLpattern,PLpattern,nrdiff)

  # rename columns
  labels <- c(
    "FP_PC1", "FP_PC2", "FN_PC1", "FN_PC2",
    "FUSION_FP", "FUSION_FN", "Tucker",
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
#first clear environment and read in MPLUS files
results <- list()
k <- 1:40
files <- sprintf("MPLUS_result%d.RData", k)
#mplus <- setNames(vector("list", length(files)), paste0("mplus", k))
for (i in seq_along(files)) {
  e <- new.env(parent = emptyenv())
  load(files[i], envir = e)          # load objects into a temporary environment
  results[[i]] <- as.list(e)           # store all loaded objects as a named list
}
RESULT <- join_simulation_results(results)$RESULT

##per condition
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
TABLE_percondition[TABLE_percondition[,6]==4,]

##following table in De Roover & Vermunt
TABLE <- matrix(nrow=16,ncol=3)  #%datasets no FP/FN for fusion and simple structure, everything correct
colnames(TABLE) <- c('FUSION', 'SIMPLE S.', 'BOTH')

#find number of datasets 100%correct: OVERALL
fusion_ISF <- sum(rowSums(RESULT[,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[,1:6]==0)==6)#ISF
TABLE[16,] <- c(fusion_ISF,simple_ISF,both_ISF)/8000

#find number of datasets 100%correct: Per nr of components
Rvec <- c(2,4)
Gvec <- c(2,4)
for (g in 1:2){
  sel <- which(RESULT$G==Gvec[g])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[g,] <- c(fusion_ISF,simple_ISF,both_ISF)/4000
}
#find number of datasets 100%correct: Per sample size
Nvec <- c(50,200,600,1000)
for (n in 1:4){
  sel <- which(RESULT$n==Nvec[n])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[2+n,] <- c(fusion_ISF,simple_ISF,both_ISF)/2000
}
#find number of datasets 100%correct: Per group size
for (q in 1:2){
  sel <- which(RESULT$R==Rvec[q])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[6+q,] <- c(fusion_ISF,simple_ISF,both_ISF)/4000
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
#find number of datasets 100%correct: Per group size
Dvec <- c(4,16)
for (d in 1:2){
  sel <- which(RESULT$NrDiff==Dvec[d])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[13+d,] <- c(fusion_ISF,simple_ISF,both_ISF)/4000
}
write.table(TABLE, file='MPLUSsimresults.txt',sep = '\t',dec = ",")

