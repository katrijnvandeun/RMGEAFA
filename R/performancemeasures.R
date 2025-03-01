performancemeasures <- function(selmodel, P1, P2, R, G) {
  #obtain best permutation based on difference in loading for P1 and P2;
  #account for sign invariance
  perm <- gtools::permutations(R, R)
  absdiff <- c()
  for (p in 1:nrow(perm)) {
    #P1
    corsign <- sign(diag(cor(P1, selmodel$loadings[[1]][, perm[p,]])))
    L1 <- (selmodel$loadings[[1]][, perm[p,]]) %*% diag(corsign)
    diffP1 <- sum(rowSums(abs(P1 - L1)))
    #P2
    corsign <- sign(diag(cor(P2, selmodel$loadings[[2]][, perm[p,]])))
    L2 <- (selmodel$loadings[[2]][, perm[p,]]) %*% diag(corsign)
    diffP2 <- sum(rowSums(abs(P2 - L2)))
    absdiff[p] <- abs(diffP1) + abs(diffP2)
  }
  bestperm <- which.min(absdiff)
  L1 <- selmodel$loadings[[1]][, perm[bestperm,]]
  L2 <- selmodel$loadings[[2]][, perm[bestperm,]]
  
  #Calculate performance measures
  tol <- 0.01 #tolerance for no difference / zero
  FP_sparse1 <-  sum(P1 == 0 & abs(round(L1, 2)) >= tol)
  FN_sparse1 <-  sum(P1 != 0 & abs(round(L1, 2)) < tol)
  FP_sparse2 <-  sum(P2 == 0 & abs(round(L2, 2)) >= tol)
  FN_sparse2 <-  sum(P2 != 0 & abs(round(L2, 2)) < tol)
  FP_fuse1 <-   sum(P1 != P2 & abs(round(L1 - L2, 2)) < tol)
  FN_fuse1 <-  sum(P1 == P2 & abs(round(L1 - L2, 2)) > tol)
  spcr_tucongrT1 <-  sum(diag(abs(psych::factor.congruence(P1, L1)))) / R
  spcr_tucongrT2 <-  sum(diag(abs(psych::factor.congruence(P2, L2)))) / R
  spcr_tucongrT <- mean(c(spcr_tucongrT1, spcr_tucongrT2))
      
  parameters <- list(
  c(
    FP_sparse1,
    FP_sparse2,
    FN_sparse1,
    FN_sparse2,
    FP_fuse1,
    FN_fuse1,
    spcr_tucongrT
  ),
  perm[bestperm, ]
)
return(parameters)
}