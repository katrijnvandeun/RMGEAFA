#' Multistart
#'
#' INSERT DESCRIPTION
#'
#' @param DATA
#' @param R
#' @param CARD
#' @param LAMBDA
#' @param SP
#' @param ALPH
#' @param MAXITER
#' @param EPS
#'
#' @return
#' @export
MULTISTART <- function(DATA, R, CARD, LAMBDA, SP, ALPH, MAXITER, EPS) {
  Loss <- Inf
  for (a in 1:length(ALPH)) {
    for (s in 1:length(SP)) {
      result <- jointSPCA(DATA, R, SP = SP[s], ALPH = ALPH[a], lambda = LAMBDA, CARD = CARD, MAXITER, EPS)
      if (result$Residual < Loss) {
        solution <- result
        Loss <- result$Residual
      }
    }
  }
  return(solution)
}
