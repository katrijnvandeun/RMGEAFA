#' SeqStrategy
#'
#' INSERT DESCRIPTION
#'
#' @param DATA
#' @param R
#' @param SP
#' @param ALPH
#' @param THR
#' @param TOL
#'
#' @return
#' @export
seqstrategy <- function(DATA, R, SP = 0.9, ALPH = 0.1, THR = 0.05, TOL = 2) {
  lambdalow <- 1e-4 # limited fusion when selecting cardinality
  J <- dim(DATA[[1]])[2]
  eps <- 10^-6

  # first -to prioritize them- check solutions with equal number of nonzero-loadings factor
  cardvec <- 3:J
  K <- length(cardvec)
  avec <- matrix(nrow = length(cardvec), ncol = 8)
  colnames(avec) <- c("lambda", "K", "IS", "PEV", "Prop0", "PropUnique", "MinNonZeroL", "MaxSDL")
  for (k in 1:K) {
    a <- IS_JSPCA(DATA, R = R, SP, ALPH, lambda = lambdalow, card = rep(cardvec[k], R), MaxIter = 20, eps)
    avec[k, ] <- c(
      lambdalow, cardvec[k], a$value, a$vaf, a$propzero, a$propunique,
      a$smallestP, a$maxsdP
    )
  }
  round(avec, 3)

  # second check solutions with cardinality at level of the loading matrix
  cardvec <- (2 * R):(J * R)
  L <- length(cardvec)
  avec <- rbind(avec, matrix(nrow = length(cardvec), ncol = 8))

  for (l in 1:L) {
    a <- IS_JSPCA(DATA, R = R, SP, ALPH, lambda = lambdalow, card = cardvec[l], MaxIter = 20, eps)
    avec[K + l, ] <- c(
      lambdalow, cardvec[l], a$value, a$vaf, a$propzero, a$propunique,
      a$smallestP, a$maxsdP
    )
  }
  round(avec, 4)

  # selection, first based on matching cardinality to unpenalized/constrained sparseness
  placeholder <- seq(1:(K + L))
  indexnonzeroL <- placeholder[avec[, "MinNonZeroL"] < THR][1] - 1 # THR holding selection
  if (indexnonzeroL < K) {
    placeholder2 <- placeholder[(K + 1):(K + L)]
    indexnonzeroL2 <- placeholder2[avec[(K + 1):(K + L), "MinNonZeroL"] < THR][1] - 1 # THR holding selection
    if (avec[indexnonzeroL2, 2] < avec[indexnonzeroL, 2] * R) {
      indexnonzeroL <- indexnonzeroL2
    } else if (avec[indexnonzeroL2, 2] == avec[indexnonzeroL, 2] * R) {
      if (round(avec[indexnonzeroL2, 4], THR) > round(avec[indexnonzeroL, 4], THR)) {
        indexnonzeroL <- indexnonzeroL2
      }
    }
  }
  selcardinality <- avec[indexnonzeroL, 2]
  if (indexnonzeroL < K + 1) {
    selcardinality <- c(rep(avec[indexnonzeroL, 2], R))
  }
  # second, to account for not finding a solution in the previous step and/or
  # occurrence of sparser solutions with almost same PVE
  samepve <- round(avec[, 4], TOL) == round(avec[indexnonzeroL, 4], TOL)
  # what follows is selection by maximizing IS among solutions with same PVE as max PVE
  # if this pve > pve of selected solution based on smallest nonzero loading
  index <- indexnonzeroL
  if (is.na(indexnonzeroL) || max(round(avec[samepve, 3], 3)) > round(avec[indexnonzeroL, 3], 3)) {
    index <- placeholder[avec[, 3] == max(avec[samepve, 3])][1]
    selcardinality <- avec[index, 2]
    if (index < K) {
      selcardinality <- rep(avec[index, 2], R)
    }
  }
  selcardinality
  AVEC <- round(avec, 4)

  lambdavec <- c(1e-4, 1e-3, 0.005, 0.0075, 0.01, 0.02, 0.025) # , 0.1, 0.2, 0.5, 1) # highest value = diff of this size is fused
  avec <- matrix(nrow = length(lambdavec), ncol = 8)
  colnames(avec) <- c("lambda", "K", "IS", "PEV", "Prop0", "PropUnique", "MinNonZeroL", "MaxSDL")
  for (k in 1:length(lambdavec)) {
    a <- IS_JSPCA(DATA, R = R, SP, ALPH, lambda = lambdavec[k], card = selcardinality, MaxIter = 20, eps)
    avec[k, ] <- c(
      lambdavec[k], selcardinality[1], a$value, a$vaf, a$propzero,
      a$propunique, a$smallestP, a$maxsdP
    )
  }
  AVEC <- rbind(AVEC, round(avec, 4))
  index <- 1:length(lambdavec)
  #index <- index[avec[, "MinNonZeroL"] >= THR]
  index <- index[which.max(avec[index, 3])]
  sellambda <- avec[index, 1]
  selmodel <- jointSPCA(DATA, R = R, SP, ALPH, lambda = sellambda, CARD = selcardinality, MaxIter = 20, eps)

  result <- list()
  result$cardinality <- selcardinality
  result$lambda <- sellambda
  result$selmodel <- selmodel
  result$details <- AVEC
  return(result)
}
