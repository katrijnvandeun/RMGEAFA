
#2 groups and 2 components
#primary loading shifts
#no of loading differences = 4
setwd("C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/cSPCA/R/Simulation")
source("../CoreFunctions/JointSPCA_finalizing.R")
source("../CoreFunctions/IS_JSPCA.R")
source("../CoreFunctions/InitialLoadings.R")
source("../CoreFunctions/seqstrategy.R")
source("../CoreFunctions/performancemeasures.R")
source("../CoreFunctions/pstr.R")

library(MASS)
n <- c(50, 200, 600, 1000)
final_result <- data.frame()
set.seed(45468)
p1g1 <- c(0,rep(sqrt(.6),9),rep(0,10))
p2g1 <- c(sqrt(.6),rep(0,9),rep(sqrt(.6),10))
P1 <- cbind(p1g1,p2g1)
p1g2 <- c(rep(sqrt(.6),2),0,rep(sqrt(.6),7),rep(0,10))
p2g2 <- c(0,0,sqrt(.6),rep(0,7),rep(sqrt(.6),10))
P2 <- cbind(p1g2,p2g2)
PSI1 <- diag(1-rowSums(P1^2))
SIGMA1 <- P1%*%t(P1)+PSI1 
PSI2 <- diag(1-rowSums(P2^2))
SIGMA2 <- P2%*%t(P2)+PSI2
for (i in 1:length(n)){
  getValue <- function(){
    data1 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA1, empirical = TRUE)
    data2 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA2, empirical = TRUE)
    DATA <- list(data1,data2)
    
    G <- 2
    R <- 2
    
    #model selection
    selmodel <- seqstrategy(DATA,R, SP=0.9, ALPH=0.1)
    
    #calculation performance measures
    parameters <- performancemeasures(selmodel,P1,P2,R,G)[[1]]    
    #round(cbind(P1,selmodel$loadings[[1]],P2,selmodel$loadings[[2]]),2)
    return(parameters)
  }
  result <- as.data.frame(t(replicate(50,getValue())))
  result$size <- c(rep(n[i],50))
  result$nrpcs <- 2
  result$nrgroups <- 2
  final_result <- rbind(final_result,result)
}

#colnames(final_result) <- c("Recovery","Tucker","Group_size")
save(final_result, file = "final_result1v2.Rda")
