#2 groups and 2 components
#crossloadings 0.4
#no of loading differences = 4
#source("C:/Users/TSB-MTO/Documents/Tra/JointSPCA_finalizing.R")
#source("C:/Users/TSB-MTO/Documents/Tra/IS_JSPCA.R")
#setwd("C:/Users/TSB-MTO/Documents/Tra/dataframes")
#setwd("C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/cSPCA/R/Simulation")
library(MASS)
n <- c(50, 200, 600, 1000)
final_result3 <- data.frame()
set.seed(45468)
p1g1 <- c(rep(sqrt(.6),10),.4,rep(0,9))
p2g1 <- c(.4,rep(0,9),rep(sqrt(.6),10))
P1 <- cbind(p1g1,p2g1)
p1g2 <- c(rep(sqrt(.6),10),rep(0,4),.4,rep(0,5))
p2g2 <- c(rep(0,4),.4,rep(0,5),rep(sqrt(.6),10))
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
  
  result3 <- as.data.frame(t(replicate(50,getValue())))
  result3$size <- c(rep(n[i],50))
  result3$nrpcs <- 2
  result3$nrgroups <- 2
  final_result3 <- rbind(final_result3,result3)
}

#colnames(final_result3) <- c("Recovery","Tucker","Group_size")
save(final_result3, file = "final_result3v2.Rda")

