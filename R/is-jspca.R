#' IS JSPCA
#'
#' IS for MODEL SELECTION regularized MGAEFA.
#' IS(K,lambda) = PVEnonregul*PVEregul*(df/(G*J*R)),
#' with    PVEnonregul the PVE in the ordinary MGAEFA
#'         PVEregul the PVE in the regularized MGAEFA with given cardinality K and fusion penalty lambda
#'         df the effective degrees of freedom, this is number of unique nonzero loadings
#'         G*J*R is the total number of loadings
#' Our modified index builds on related work for the BIC for the fused lasso.
#' See papers of Guo et al. (2010) on sparse fused lasso for more references; also check
#' Nowak et al. (sparse fused latent features, 2011)
#'
#' Note: Original Index of sparseness = VAFunregul*VAFregul*(df/JR)
#' Originally proposed by Gajjar et al. 2017
#'
#' @param X
#' @param R
#' @param SP
#' @param ALPH
#' @param lambda
#' @param card
#' @param MaxIter
#' @param eps
#'
#' @return
#' @export
IS_JSPCA <- function(X, R, SP, ALPH, lambda, card, MaxIter, eps) {
  # setting default values for some parameters
  if (is.null(SP)) SP <- 0.9
  if (is.null(ALPH)) ALPH <- 0.1
  TOL <- 3
  if (is.null(MaxIter)) MaxIter <- 20
  if (is.null(eps)) eps <- 1e-4
  # initialize list for output
  IS <- list()
  G <- length(X)
  J <- dim(X[[1]])[2]
  # 1. Obtain Joint SPCA solution
  JSPCA <- jointSPCA(X, R, SP, ALPH, lambda = lambda, CARD = card, MaxIter, eps)
  # 2. Calculate PVE and create G x (JR) loading matrix
  loadingMAT <- matrix(data = NA, nrow = G, ncol = J * R)
  sspca <- 0
  sstotal <- 0
  for (g in 1:G) {
    sstotal <- sstotal + sum(rowSums(X[[g]]^2))
    d <- svd(X[[g]])$d
    sspca <- sspca + sum((d[1:R])^2)
    loadingMAT[g, ] <- c(JSPCA$loadings[[g]])
  }
  ssjspca <- sstotal - JSPCA$Residual
  # 3. Calculation degrees of freedom
  dfjspca <- 0
  # obtain nr of unique nonzero loadings
  loadingMAT_rounded <- round(loadingMAT, TOL) # rounding needed else all diff
  for (jr in 1:(J * R)) {
    dfjspca <- dfjspca + length(unique(loadingMAT_rounded[, jr][loadingMAT[, jr] != 0]))
  }
  nrcoef <- J * R * G
  nrzero <- nrcoef - (G * sum(card))
  nrzeqcoef <- nrcoef - dfjspca
  IS$value <- sspca * ssjspca / sstotal^2 * nrzeqcoef / nrcoef
  IS$vaf <- 1 - (JSPCA$Residual / sstotal)
  IS$propunique <- dfjspca / (nrcoef - nrzero)
  IS$propzero <- nrzero / nrcoef
  IS$smallestP <- ifelse(sum(rowSums(loadingMAT != 0)) < G * sum(card),
    0, min(abs(loadingMAT[loadingMAT != 0]))
  )
  IS$maxsdP <- max(apply(loadingMAT, 2, sd))
  return(IS)
}
