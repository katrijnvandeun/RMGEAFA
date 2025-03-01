#### Title: Joint Sparse PCA
#### Author: Katrijn Van Deun, PhD.
#### Created: September 6, 2021
#### Last modified: February 4, 2024

#########################################################################################
#####                                 Preliminaries                                 #####
#########################################################################################
# Let Xg be a Ig x J matrix containing the data for group g, with g = 1, ..., G; note
# that for all groups the same set of J variables is assumed
# Joint sparse PCA solves the following optimization problem:
#     min_{Tg,Pg} \sum_g ||X_g-T_gP_g^T||^2 + \lambda\sum_{j,r}\sum_{g,g'}|p_{j,r}^{(g)}-p_{j,r}^{(g')}|)
# subject to Tg^TTg = I and a cardinality constraint on the loadings (either at 
# the level of the matrix or per component).

# Clear workspace
#rm(list = ls())

# Packages
library(matlib)          # inverse of a matrix inv()



#########################################################################################
#####                               JOINT Sparse PCA                                #####
#########################################################################################

#Input: DATA, a list containing the data matrix of each group,
#R the number of components, lambda (between 0 and 2) the strength of the penalties,
#CARD the number of nonzero loadings, and MaxIter the maximum number of iterations before stopping

jointSPCA <- function(DATA, R, SP, ALPH, lambda, CARD, MaxIter, eps) {
  #CASE 1: SEPARATE PCAs WHEN CARD == 0 & LAMBDA == 0
  if (lambda == 0 & sum(CARD) == 0){ 
    G <- length(DATA)
    ssX <- 0
    residual <- 0
    for (g in 1:G){
      svddata <- svd(DATA[[g]])
      ssX <- ssX+sum(rowSums(DATA[[g]]^2))
      loadings <- svddata$v[,1:R]%*%diag(svddata$d[1:R])/sqrt(dim(DATA[[g]])[1])
      scores <- svddata$u[, 1:R]*sqrt(dim(DATA[[g]])[1])
      XHAT <- scores%*%t(loadings)
      residual <- residual+sum(rowSums((XHAT-DATA[[g]])^2))
    }
    Loss <- residual/ssX
    result <- list()
    result$loadings <- loadings
    result$scores <- scores
    result$Loss <- Loss
  } else if (lambda == 0 & sum(CARD) > 0) {
    #CASE 2: USLPCA WHEN CARD > 0 & LAMBDA == 0
    results <- USLPCA(DATA, R, TYPE, CARD, MaxIter, eps) 
    #loadings <- results$loadings
    #scores <- results$scores
    #Lossvec <- results$Lossvec
  } else {
    #CASE 3: JOINT SPCA WHEN CARD > 0 & LAMBDA > 0
    results <- JSPCA(DATA, R, SP, ALPH, lambda, CARD, MaxIter, eps)
    #loadings <- results$loadings
    #scores <- results$scores
    #Lossvec <- results$Lossvec
  }
  return(results)
}  
  

###############################################################################
############################SUBFUNCTIONS#######################################
###############################################################################

##1. LOSS
RESIDUAL <- function(DATA, SCORES, LOADINGS){
  G <- length(DATA)
  res <- 0
  for (g in 1:G) {
    XHAT <- SCORES[[g]]%*%t(LOADINGS[[g]])
    res <- res+sum(rowSums((XHAT-DATA[[g]])^2))
  }
  return(res)
}


##2. USLPCA
USLPCA <- function(DATA, R, SP, ALPH, CARD, MaxIter, eps) {
  G <- length(DATA)
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc        #Used to compute convergence criterium
  loadings <- INITIALLOADINGS(DATA, R, CARD, SP, ALPH)
  A <- loadings
  scores <- matrix(nrow = Ivec[1], ncol = R)
  scores <- lapply(1:G, function(x) scores) 
  if (sum(length(CARD)) > 1){
    CARDc <- J-CARD #Number of loadings -per PC- that have to be zero
    while (convAO == 0) {
      #1. Update component scores
      for (g in 1:G) {
        XP <- DATA[[g]] %*% loadings[[g]] / (sqrt(Ivec[g]))
        scores[[g]] <- sqrt(Ivec[g]) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
        A[[g]] <- t(DATA[[g]]) %*% scores[[g]] / Ivec[g]
      }
      #Calculate loss
      Loss <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Loss)
      #2. Update loadings
      for (g in 1:G) {
        loadings[[g]] <- A[[g]]
        for (r in 1:R){
          ind <- sort(abs(A[[g]][,r]), index.return = TRUE)
          loadings[[g]][ind$ix[1:CARDc[r]],r] <- 0
        } 
      }
      #Calculate loss
      Lossu <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Lossu)
      if (iter > MaxIter) {
        convAO <- 1
      }
      if (abs(Lossc-Lossu) < eps) {
        convAO <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }
  } else {
    CARDc <- J*R-CARD #Number of loadings -over all PCs- that have to be zero
    while (convAO == 0) {
      #1. Update component scores
      for (g in 1:G) {
        XP <- DATA[[g]] %*% loadings[[g]] / (sqrt(Ivec[g]))
        scores[[g]] <- sqrt(Ivec[g]) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
        A[[g]] <- t(DATA[[g]]) %*% scores[[g]] / Ivec[g]
      }
      #Calculate loss
      Loss <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Loss)
      #2. Update loadings
      for (g in 1:G) {
        loadings[[g]] <- A[[g]]
        ind <- sort(abs(A[[g]]), index.return = TRUE)
        loadings[[g]][ind$ix[1:CARDc]] <- 0
      }
      #Calculate loss
      Lossu <- RESIDUAL(DATA,scores,loadings)
      Lossvec <- c(Lossvec,Lossu)
      if (iter > MaxIter ) {
        convAO <- 1
      }
      if (abs(Lossc-Lossu) < eps) {
        convAO <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }
  }
  uslpca <- list('scores' = scores, 'loadings' = loadings, 'Lossvec' = Lossvec, 'Residual' = Loss)
  return(uslpca)
}

##3. JSPCA
JSPCA <- function(DATA, R, SP, ALPH, lambda, CARD, MaxIter, eps){
  #default values for ADMM TUNING PARAMETERS
  rho = 1
  nu = 1
  ##Constants associated to the data
  G <- length(DATA)
  J <- dim(DATA[[1]])[2]
  Ivec <- vector(mode = 'numeric', length = G)
  Ivec[1] <- dim(DATA[[1]])[1]
  CONCdata <- DATA[[1]]
  scores <- list()
  scores[[1]] <- matrix(nrow = Ivec[1], ncol = R)
  if(G>1){
  for (g in 2:G) {
    CONCdata <- rbind(CONCdata, DATA[[g]])
    Ivec[g] <- dim(DATA[[g]])[1]
    scores[[g]] <- matrix(nrow = Ivec[g], ncol = R)
  }
  }
  ssX <- sum(rowSums(CONCdata ^ 2))
  ##Constants in the iterations
  D <- DIFFOP(G)
  invD <- inv(diag(G) + rho * t(D) %*% D)#
  ##INITIALIZATION LOADINGS
  #loadings <- INITLOADINGS(DATA, R, CARD, TYPE)
  loadings <- INITIALLOADINGS(DATA, R, CARD, SP, ALPH)
  ##PRE-ALLOCATION MEMORY FOR INNER ADMM: GxJR matrices
  QMAT <- matrix(nrow = G, ncol = J * R)
  for (g in 1:G){
    QMAT[g,] <- c(loadings[[g]])
  }
  pen <- lambda * (sum(rowSums(abs(D %*% QMAT))))
  YMAT <- QMAT
  AMAT <- QMAT
  PMAT <- QMAT
  A <- loadings
  #LAGRANGE MULTIPLIERS INNER LOOP
  WMAT <- matrix(1, nrow = G * (G - 1) / 2, ncol = J * R, byrow = TRUE)
  #Initialize U [LAGRANGE MULTIPLIERS OUTER LOOP] (JxR)
  onemat <- matrix(1, nrow = J, ncol = R)
  U <- lapply(1:G, function(r)  onemat)#outer
  innervec <- c()         #keep track of number of iterations inner admm
  outervec <- c()         #keep track of number of iterations outer admm
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc
  if (sum(length(CARD)) > 1){
    CARDc <- J-CARD
    while (convAO == 0) {
      ############UPDATE SCORE MATRICES Tg####################
      for (g in 1:G) {
        XP <- DATA[[g]] %*% loadings[[g]] / (sqrt(Ivec[g]))
        scores[[g]] <- sqrt(Ivec[g]) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
        A[[g]] <- t(DATA[[g]]) %*% scores[[g]] / Ivec[g]
        AMAT[g, ] <- c(A[[g]])
      }
      residual <- RESIDUAL(DATA,scores,loadings)
      pen <- lambda * (sum(rowSums(abs(D %*% PMAT))))
      Lossu <- (0.5 * residual + pen) / (2*ssX)
      Lossvec <- c(Lossvec,Lossu)
      ############UPDATE LOADING MATRICES Pg####################
      ###OUTER ADMM###
      convadmmouter <- 0 #convergence outer admm loop
      iterouter <- 0
      while (convadmmouter == 0 ) {
        QMATold <- QMAT #needed to check convergence outer ADMM
        #STEP1 OUTER ADMM: Update primal P
        for (g in 1:G) {
          part2 <- nu * matrix(PMAT[g, ], nrow = J, ncol = R) + nu * U[[g]]
          loadings[[g]] <- (A[[g]] + part2) / (1 + nu)
          Pest <-  (A[[g]]) ^ 2 + (part2) ^ 2#sorts from smallest to largest
          for (r in 1:R){
            ind <- sort(Pest[,r], index.return = TRUE)
            loadings[[g]][ind$ix[1:CARDc[r]],r] <- 0
          }
          PMAT[g, ] <- c(loadings[[g]])
          YMAT[g, ] <- c(loadings[[g]] - U[[g]])
        }
        #STEP2 OUTER ADMM: Update dual Q
        ##INNER ADMM = FUSED LASSO###
        ZMAT = D %*% PMAT
        convadmminner <- 0
        iterinner <- 0
        while (convadmminner == 0) {
          #1. Update primal Q
          QMAT <- invD %*% (YMAT + (rho * t(D) %*% (ZMAT - WMAT)))
          DIFF <- D %*% QMAT
          a <- DIFF + WMAT
          #2. Update dual Z
          ZMATup <- sign(a) * mapply(max, abs(a) - (lambda / rho), 0)
          # 3.Update lagrange multipliers W
          WMAT <- WMAT + DIFF - ZMAT
          primres <- sum(rowSums((DIFF-ZMATup)^2))/sum(rowSums(QMAT^2))
          dualres <- sum(rowSums((ZMATup-ZMAT)^2))/sum(rowSums(WMAT^2))
          if (primres < eps & dualres < eps) {
            convadmminner <- 1
          }
          iterinner <- iterinner + 1
          if (iterinner > MaxIter) {
            convadmminner <- 1
          }
          ZMAT <- ZMATup
          iterinner <- iterinner+1
        }
        innervec <- c(innervec,iterinner)
        #STEP3 OUTER ADMM: Update Lagrange multipliers U
        ssqU <- 0
        for (g in 1:G) {
          U[[g]] <- U[[g]] + nu * (matrix(QMAT[g, ], nrow = J, ncol = R) - loadings[[g]])
          ssqU <- ssqU + sum(rowSums(U[[g]] ^ 2))
        }
        primres <- sum(rowSums((QMAT - PMAT) ^ 2)) / sum(rowSums(QMAT^2))
        dualres <- sum(rowSums((QMAT - QMATold) ^ 2)) / ssqU
        if (primres < eps & dualres < eps) {
          convadmmouter <- 1
        }
        if (iterouter > MaxIter) {
          convadmmouter <- 1
        } 
        iterouter <- iterouter + 1
      }
      outervec <- c(outervec,iterouter)
      residual <- RESIDUAL(DATA,scores,loadings)
      pen <- lambda * (sum(rowSums(abs(D %*% PMAT))))
      Lossu <- (0.5 * residual + pen) / (2*ssX)
      Lossvec <- c(Lossvec, Lossu)
      if (iter > MaxIter ) {
        convAO <- 1
      }
      if (abs(Lossc-Lossu) < eps) {
        convAO <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }
  } else {
    CARDc <- J*R-CARD
    while (convAO == 0) {
      ############UPDATE SCORE MATRICES Tg####################
      for (g in 1:G) {
        XP <- DATA[[g]] %*% loadings[[g]] / (sqrt(Ivec[g]))
        scores[[g]] <- sqrt(Ivec[g]) * svd(XP, R, R)$u %*% t(svd(XP, R, R)$v)
        A[[g]] <- t(DATA[[g]]) %*% scores[[g]] / Ivec[g]
        AMAT[g, ] <- c(A[[g]])
      }
      residual <- RESIDUAL(DATA,scores,loadings)
      pen <- lambda * (sum(rowSums(abs(D %*% PMAT))))
      Lossu <- (0.5 * residual + pen) / (2*ssX)
      Lossvec <- c(Lossvec, Lossu)
      ############UPDATE LOADING MATRICES Pg####################
      ###OUTER ADMM###
      convadmmouter <- 0 #convergence outer admm loop
      iterouter <- 0
      while (convadmmouter == 0 ) {
        QMATold <- QMAT #needed to check convergence outer ADMM
        #STEP1 OUTER ADMM: Update primal P
        for (g in 1:G) {
          part2 <- nu * matrix(PMAT[g, ], nrow = J, ncol = R) + nu * U[[g]]
          loadings[[g]] <- (A[[g]] + part2) / (1 + nu)
          Pest <-  (A[[g]]) ^ 2 + (part2) ^ 2#sorts from smallest to largest
          ind <- sort(Pest, index.return = TRUE)
          loadings[[g]][ind$ix[1:CARDc]] <- 0
          PMAT[g, ] <- c(loadings[[g]])
          YMAT[g, ] <- c(loadings[[g]] - U[[g]])
        }
        #STEP2 OUTER ADMM: Update dual Q
        ##INNER ADMM = FUSED LASSO###
        ZMAT = D %*% PMAT
        convadmminner <- 0
        iterinner <- 0
        while (convadmminner == 0) {
          #1. Update primal Q
          QMAT <- invD %*% (YMAT + (rho * t(D) %*% (ZMAT - WMAT)))
          DIFF <- D %*% QMAT
          a <- DIFF + WMAT
          #2. Update dual Z
          ZMATup <- sign(a) * mapply(max, abs(a) - (lambda / rho), 0)
          # 3.Update lagrange multipliers W
          WMAT <- WMAT + DIFF - ZMAT
          primres <- sum(rowSums((DIFF-ZMATup)^2))/sum(rowSums(QMAT^2))
          dualres <- sum(rowSums((ZMATup-ZMAT)^2))/sum(rowSums(WMAT^2))
          if (primres < eps & dualres < eps) {
            convadmminner <- 1
          }
          iterinner <- iterinner + 1
          if (iterinner > MaxIter) {
            convadmminner <- 1
          }
          ZMAT <- ZMATup
          iterinner <- iterinner+1
        }
        innervec <- c(innervec,iterinner)
        #STEP3 OUTER ADMM: Update Lagrange multipliers U
        ssqU <- 0
        for (g in 1:G) {
          U[[g]] <- U[[g]] + nu * (matrix(QMAT[g, ], nrow = J, ncol = R) - loadings[[g]])
          ssqU <- ssqU + sum(rowSums(U[[g]] ^ 2))
        }
        primres <- sum(rowSums((QMAT - PMAT) ^ 2)) / sum(rowSums(QMAT^2))
        dualres <- sum(rowSums((QMAT - QMATold) ^ 2)) / ssqU
        if (primres < eps & dualres < eps) {
          convadmmouter <- 1
        }
        if (iterouter > MaxIter) {
          convadmmouter <- 1
        } 
        iterouter <- iterouter + 1
      }
      outervec <- c(outervec,iterouter)
      residual <- RESIDUAL(DATA,scores,loadings)
      pen <- lambda * (sum(rowSums(abs(D %*% PMAT))))
      Lossu <- (0.5 * residual + pen) / (2*ssX)
      Lossvec <- c(Lossvec, Lossu)
      if (iter > MaxIter ) {
        convAO <- 1
      }
      if (abs(Lossc-Lossu) < eps) {
        convAO <- 1
      }
      iter <- iter + 1
      Lossc <- Lossu
    }  
  }
  return_jspca <- list()
  return_jspca$loadings <- loadings
  return_jspca$scores <- scores
  return_jspca$Lossvec <- Lossvec
  #return_varselect$Residual <- residual
  return_jspca$iterinner <- innervec
  return_jspca$iterouter <- outervec
  return_jspca$Residual <- residual
  return(return_jspca)
}



##4. INVERSE DIFFERENCE OPERATOR
DIFFOP <- function(G){
  ##difference operator
  D <- matrix(0, nrow = G * (G - 1) / 2, ncol = G)
  counter <- 1
  for (i in 1:(G - 1)) {
    for (j in (i + 1):G) {
      D[counter, i] <- 1
      D[counter, j] <- -1
      counter <- counter + 1
    }
  }
  return(D)
}
##5. INITIAL LOADINGS: FIRST VERSION, PRINCIPAL AXES BASED
#TYPE = 'random' is sca with random zero positions defined per group
#       'pca-sep' is separate pcas with rotation to congruence to SCA solution followed by thresholding
#       'sca' is same loading matrix for all groups, zeros by thresholding
INITLOADINGS <- function(DATA, R, CARD, TYPE){
  G <- length(DATA)
  J <- dim(DATA[[1]])[2]
  CONCdata <- DATA[[1]]
  Ivec <- vector(mode = 'numeric', length = G)
  Ivec[1] <- dim(DATA[[1]])[1]
  for (g in 2:G) {
    CONCdata <- rbind(CONCdata, DATA[[g]])
    Ivec[g] <- dim(DATA[[g]])[1]
  }
  ##initialize list of G loading matrices
  #low cardinality => spread zeros over components if cardinality is matrixwise
  if (length(CARD)==1){
    CARD=rep(floor(CARD/R),R)
  }
  #1. 'random'
  if (TYPE == 'random'){
    svdCONC <- svd(CONCdata, R, R)
    loadings <- svdCONC$v[,1:R] %*% diag(svdCONC$d[1:R]) / sqrt(Ivec[1])#let loadings reflect correlation of variables with Tconc
    loadings <- lapply(1:G, function(x)  loadings)
    if (length(CARD)>1){
      CARDc <- J-CARD
      for (g in 1:G) {
        for (r in 1:R){
          rands <- sample(J, CARDc[r])#without replacement
          loadings[[g]][rands,r] <- 0
          loadings[[g]] <- loadings[[g]]*sqrt(Ivec[1])/sqrt(Ivec[g])
        } 
      }
    } else {
      CARDc <- J*R-CARD
      for (g in 1:G){
        rands <- sample(J*R, CARDc)#without replacement
        loadings[[g]][rands] <- 0
        loadings[[g]] <- loadings[[g]]*sqrt(Ivec[1])/sqrt(Ivec[g])
      }
    }
  } else if (TYPE == 'sca') {
    svdCONC <- svd(CONCdata, R, R)
    loadingsconc <- svdCONC$v[,1:R] %*% diag(svdCONC$d[1:R]) / (sqrt(G) * sqrt(sum(Ivec)))#let loadings reflect correlation of variables with Tconc
    P <- loadingsconc
    if (length(CARD)>1){
      CARDc <- J-CARD
      for (r in 1:R){
        ind <- sort(abs(P[,r]), index.return = TRUE)
        P[ind$ix[1:CARDc[r]],r] <- 0
      } 
    } else {
      CARDc <- J*R-CARD
      ind <- sort(abs(P), index.return = TRUE)
      P[ind$ix[1:CARDc]] <- 0
    }
    loadings <- lapply(1:G, function(x)  P)
    for (g in 1:G){
      loadings[[g]] <- loadings[[g]]*sqrt(Ivec[1])/sqrt(Ivec[g])
    }
  } else { #'pca-sep'
    svdCONC <- svd(CONCdata, R, R)
    loadingsconc <- svdCONC$v[,1:R] %*% diag(svdCONC$d[1:R]) / (sqrt(G) * sqrt(sum(Ivec)))#let loadings reflect correlation of variables with Tconc
    scores <- lapply(DATA, function(x,R) matrix(nrow = dim(x)[1],ncol=R), R)
    svdresults <- lapply(DATA, svd)
    loadings <- lapply(svdresults,function(x,R) x$v[,1:R]%*%diag(x$d[1:R]),R)
    #scores <- lapply(svdresults,function(x,R) x$u[,1:R],R)
    for (g in 1:G) {
      loadings[[g]] <- svdresults[[g]]$v[,1:R] %*% diag(svdresults[[g]]$d[1:R]) / sqrt(Ivec[g])
      tuck <- TuckerCoef(loadingsconc, loadings[[g]])
      loadings[[g]] <- loadings[[g]][,tuck$perm]
      corsign <- sign(diag(cor(loadingsconc,loadings[[g]])))
      loadings[[g]] <- loadings[[g]]%*%diag(corsign)
      if (length(CARD)>1){
        CARDc <- J-CARD
        for (r in 1:R){
          ind <- sort(abs(loadings[[g]][,r]), index.return = TRUE)#without replacement
          loadings[[g]][ind$ix[1:CARDc[r]],r] <- 0
        } 
      } else {
        CARDc <- J*R-CARD
        ind <- sort(abs(loadings[[g]]),index.return = TRUE)#without replacement
        loadings[[g]][ind$ix[1:CARDc]] <- 0
      }
    }
  }
  ###account for permutation, reflection
  for (g in 2:G) {
    perm <- gtools::permutations(R, R)
#    corrate <- vector(length = nrow(perm))
#    for (i in 1:nrow(perm)) {
#      corrate[i]<-num_correct(loadings[[1]], loadings[[g]][,perm[i,]])
#    }
#    for (r in 1:R) {
#      corsign <- sign(sum(sign(loadings[[1]][,r])*sign(loadings[[g]][,r])))
#      loadings[[g]] <- loadings[[g]]*corsign
#    }
    svdCONC <- svd(CONCdata, R, R)
    loadingsconc <- svdCONC$v[,1:R] %*% diag(svdCONC$d[1:R]) / (sqrt(G) * sqrt(sum(Ivec)))#let loadings reflect correlation of variables with Tconc
    absdiff <- c()
    for (p in 1:nrow(perm)) {
      corsign <- sign(diag(cor(loadingsconc, loadings[[g]][, perm[p,]])))
      L1 <- (loadings[[g]][, perm[p,]]) %*% diag(corsign)
      diff <- sum(rowSums(abs(loadingsconc - L1)))
      absdiff[p] <- abs(diff)
    }
    bestperm <- which.min(absdiff)
    corsign <- sign(diag(cor(loadingsconc, loadings[[g]][, perm[bestperm,]])))
    loadings[[g]] <- (loadings[[g]][,perm[bestperm,]]) %*% diag(corsign)
    #corsign <- sign(diag(cor(loadings[[1]],loadings[[g]])))
    #loadings[[g]] <- loadings[[g]]%*%diag(corsign)
  }
  return(loadings)
}

###function Zhengguo for recovery rate
num_correct <- function (TargetP, EstimatedP){
  total_vnumber <- dim(TargetP)[1] * dim(TargetP)[2]
  TargetP[which(TargetP != 0)] <- 1
  sum_select <- sum(TargetP)
  sum_zero <- total_vnumber - sum_select
  EstimatedP[which(EstimatedP != 0)] <- 1
  total_correct <- sum(TargetP == EstimatedP) # this is the total number of variables correctedly selected and zeros correctly retained
  prop_correct <- total_correct/total_vnumber
  return(prop_correct)
}