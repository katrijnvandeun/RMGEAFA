#2 groups and 2 components
#primary loading decreases 0.4
#no of loading differences = 4
#source("C:/Users/TSB-MTO/Documents/Tra/JointSPCA_finalizing.R")
#source("C:/Users/TSB-MTO/Documents/Tra/IS_JSPCA.R")
library(MASS)
n <- c(50,200,600,1000)
final_result7 <- data.frame()
set.seed(45468)
p1g1 <- c(sqrt(.6)-.4,rep(sqrt(.6),9),rep(0,10))
p2g1 <- c(rep(0,10),sqrt(.6)-.4,rep(sqrt(.6),9))
P1 <- cbind(p1g1,p2g1)
p1g2 <- c(rep(sqrt(.6),4),sqrt(.6)-.4,rep(sqrt(.6),5),rep(0,10))
p2g2 <- c(rep(0,10),rep(sqrt(.6),4), sqrt(.6)-.4, rep(sqrt(.6),5))
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
  
  result7 <- as.data.frame(t(replicate(50,getValue())))
  result7$size <- c(rep(n[i],50))
  result7$nrpcs <- 2
  result7$nrgroups <- 2
  final_result7 <- rbind(final_result7,result7)
}

#colnames(final_result7) <- c("Recovery","Tucker","Group_size")
save(final_result7, file = "final_result7v2.Rda")