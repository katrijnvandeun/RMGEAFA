### Imports
library(MASS)

if (sys.nframe() == 0L) { # Check if the file is being sourced or run on its own for the impoert pahts to work correctly
  source("./R/joint-spca-finalizing.R")
  source("./R/seqstrategy.R")
  source("./R/multistart.R")
  source("./R/performance-measures.R")
  source("./R/pstr.R") # 2nd version: Rotation towards simple structure
  source("./R/initial-loadings.R")
  source("./R/is-jspca.R")
} else {
  source("../R/joint-spca-finalizing.R")
  source("../R/seqstrategy.R")
  source("../R/multistart.R")
  source("../R/performance-measures.R")
  source("../R/pstr.R") # 2nd version: Rotation towards simple structure
  source("../R/initial-loadings.R")
  source("../R/is-jspca.R")
}

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

simulation <- function(groups, components, differences, CL, PLdecrease, seed = 45468) {
  N <- c(50, 200, 600, 1000)
  final_result <- data.frame()
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
        DATA <- list(data1, data2)
      } else if (groups == 4) {
        data1 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data2 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
        data3 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA1, empirical = TRUE)
        data4 <- mvrnorm(n = N[i], mu = rep(0, 20), Sigma = SIGMA2, empirical = TRUE)
        DATA <- list(data1, data2, data3, data4)
      }

      G <- groups
      R <- components

      # model selection
      selmodel <- seqstrategy(DATA, R, SP = 0.9, ALPH = 0.1)
      # calculation performance measures
      parameters <- performancemeasures(selmodel$selmodel, P1, P2, R, G)[[1]]
      return(parameters)
    }
    result <- as.data.frame(t(replicate(50, getValue())))
    result$size <- c(rep(N[i], 50))
    result$nrpcs <- components
    result$nrgroups <- groups
    final_result <- rbind(final_result, result)
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
    results[[i]] <- do.call(simulation, as.list(as.numeric(comb)))
  }
  return(results)
}

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
  dataG2R2 <- rbind(results[[G2R2_ids]])
  # one big dataframe for 2 groups 4 components
  dataG2R4 <- rbind(results[[G2R4_ids]])
  # one big dataframe for 4 groups 2 components
  dataG4R2 <- rbind(results[[G4R2_ids]])
  # one big dataframe for 4 groups 4 components
  dataG4R4 <- rbind(results[[G4R4_ids]])

  # join
  result <- rbind(dataG2R2, dataG2R4, dataG4R2, dataG4R4)

  # rename columns
  labels <- c(
    "FP_PC1", "FP_PC2", "FN_PC1", "FN_PC2",
    "FUSION_FP", "FUSION_FN", "Tucker",
    "n", "R", "G", "Condition", "Diff"
  )
  colnames(result) <- labels

  # Format output
  output <- list()
  output$G2R2 <- dataG2R2
  output$G2R4 <- dataG2R4
  output$G4R2 <- dataG4R2
  output$G4R2 <- dataG4R2
  output$RESULT <- result

  return(output)
}
