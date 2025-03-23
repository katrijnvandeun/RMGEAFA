#' Initial Loadings
#'
#' Initial loading matrices are generated that are semi-rational with two
#' options for the rational configuration: either the same for all groups
#' based on the cardinality constrained SCA solution or specific for each
#' group based on the cardinality constrained PCA-SEP solutions.
#' The initial SVD based configuration is rotated towards simple structure.
#'
#' @param DATA
#' @param R
#' @param CARD
#' @param SP Signal proportion for semi-rational configurations
#' @param ALPH Weight in linear combination of SCA and PCA-SEP solutions.
#' when ALPH=0 the (perturbed) SCA solution results
#' when ALPH=1 the (perturbed) PCA-SEP solution results
#'
#' @return
#' @export
INITIALLOADINGS <- function(DATA, R, CARD, SP, ALPH) {
  G <- length(DATA)
  J <- dim(DATA[[1]])[2]
  CONCdata <- c() # DATA[[1]]
  Ivec <- vector(mode = "numeric", length = G)
  if (J > 1000) {
    # large J => eye(J) and kmeans do not work; use SCA as initial config
    Ivec[1] <- dim(DATA[[1]])[1]
    CONCdata <- DATA[[1]]
    for (g in 2:G) {
      CONCdata <- rbind(CONCdata, DATA[[g]])
      Ivec[g] <- dim(DATA[[g]])[1]
    }
    svdCONC <- svd(CONCdata, R, R)
    loadingsconc <-
      svdCONC$v[, 1:R] %*% diag(svdCONC$d[1:R]) / (sqrt(G) * sqrt(sum(Ivec))) # let loadings reflect correlation of variables with Tconc
    P <- loadingsconc
    if (length(CARD) > 1) {
      CARDc <- J - CARD
      for (r in 1:R) {
        ind <- sort(abs(P[, r]), index.return = TRUE)
        P[ind$ix[1:CARDc[r]], r] <- 0
      }
    } else {
      CARDc <- J * R - CARD
      ind <- sort(abs(P), index.return = TRUE)
      P[ind$ix[1:CARDc]] <- 0
    }
    loadings <- lapply(1:G, function(x) {
      P
    })
    for (g in 1:G) {
      loadings[[g]] <- loadings[[g]] * sqrt(Ivec[1]) / sqrt(Ivec[g])
    }
    loadings_pcasca <- loadings
  } else { # note that mvrnorm has problems in HD setting if empir=T
    for (g in 1:G) {
      Ig <- dim(DATA[[g]])[1]
      NOISE <- mvrnorm(
        n = Ig,
        mu = rep(0, J),
        Sigma = diag(rep(1, J)),
        empirical = FALSE
      )
      ssxg <- sum(rowSums(DATA[[g]]^2))
      ssnoise <- sum(rowSums(NOISE^2))
      f <- sqrt(((1 - SP) / SP) * ssxg / ssnoise)
      DATA[[g]] <- DATA[[g]] + f * NOISE
      CONCdata <- rbind(CONCdata, DATA[[g]])
      Ivec[g] <- Ig
    }
    svdCONC <- svd(CONCdata, R, R)
    loadings <-
      svdCONC$v[, 1:R] %*% diag(svdCONC$d[1:R]) / (sqrt(sum(Ivec))) # let loadings reflect correlation of variables with Tconc
    # rotate towards simple structure using partially specified target rotation
    # partition variables into K clusters based on cross-product similarity
    c <- kmeans(loadings, R, nstart = 5)
    TARGET <- matrix(0, nrow = J, ncol = R)
    for (r in 1:R) {
      TARGET[c$cluster == r, r] <- 1
    }
    WEIGHTS <- TARGET
    B <- pstr(loadings, TARGET, WEIGHTS, 50, 1e-4)
    loadings <- loadings %*% B$Bmatrix
    loadings_sca <- lapply(1:G, function(x) {
      loadings
    })
    loadings_pcasep <- list()
    loadings_pcasca <- list()
    perm <- gtools::permutations(R, R)
    svdresults <- lapply(DATA, svd)
    for (g in 1:G) {
      loadings_pcasep[[g]] <-
        svdresults[[g]]$v[, 1:R] %*% diag(svdresults[[g]]$d[1:R]) / (sqrt(Ivec[g]))
      absdiff <- c()
      for (p in 1:nrow(perm)) {
        corsign <-
          sign(diag(cor(loadings, loadings_pcasep[[g]][, perm[p, ]])))
        L1 <- (loadings_pcasep[[g]][, perm[p, ]]) %*% diag(corsign)
        diff <- sum(rowSums(abs(loadings - L1)))
        absdiff[p] <- abs(diff)
      }
      bestperm <- which.min(absdiff)
      corsign <-
        sign(diag(cor(loadings, loadings_pcasep[[g]][, perm[bestperm, ]])))
      loadings_pcasep[[g]] <-
        (loadings_pcasep[[g]][, perm[bestperm, ]]) %*% diag(corsign)
      loadings_pcasca[[g]] <-
        ALPH * loadings_pcasep[[g]] + (1 - ALPH) * loadings
    }
    # inspection
    round(cbind(loadings, loadings_pcasep[[g]], loadings_pcasca[[g]]), 2)
    # low cardinality => spread zeros over components if cardinality is matrixwise
    for (g in 1:G) {
      if (length(CARD) == 1) {
        CARD <- rep(floor(CARD / R), R)
      }
      if (length(CARD) > 1) {
        CARDc <- J - CARD
        for (r in 1:R) {
          ind <-
            sort(abs(loadings_pcasca[[g]][, r]), index.return = TRUE) # without replacement
          loadings_pcasca[[g]][ind$ix[1:CARDc[r]], r] <- 0
        }
      } else {
        CARDc <- J * R - CARD
        ind <-
          sort(abs(loadings_pcasca[[g]]), index.return = TRUE) # without replacement
        loadings_pcasca[[g]][ind$ix[1:CARDc]] <- 0
      }
    }
  }
  return(loadings_pcasca)
}
