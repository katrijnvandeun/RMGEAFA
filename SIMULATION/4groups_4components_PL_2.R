#4 groups and 4 components
#primary loading decrease 0.2
#no of loading differences = 4
#setwd("C:/Users/ASUS/surfdrive/Traineeship 1 - Tra Le/Final result dataframes")
#source("C:/Users/TSB-MTO/Documents/Tra/JointSPCA_finalizing.R")
#source("C:/Users/TSB-MTO/Documents/Tra/IS_JSPCA.R")
library(MASS)
n <- c(50,200,600,1000)
final_result38 <- data.frame()
set.seed(8790)
p1g1 <- c(sqrt(.6)-.2,rep(sqrt(.6),4), rep(0,15))
p2g1 <- c(rep(0,5),sqrt(.6)-.2,rep(sqrt(.6),4),rep(0,10))
p3g1 <- c(rep(0,10),rep(sqrt(.6),5),rep(0,5))
p4g1 <- c(rep(0,15),rep(sqrt(.6),5))
P1 <- cbind(p1g1,p2g1,p3g1,p4g1)
p1g2 <- c(rep(sqrt(.6),2),sqrt(.6)-.2,rep(sqrt(.6),2),rep(0,15))
p2g2 <- c(rep(0,5),rep(sqrt(.6),2),sqrt(.6)-.2,rep(sqrt(.6),2),rep(0,10))
p3g2 <- c(rep(0,10),rep(sqrt(.6),5),rep(0,5))
p4g2 <- c(rep(0,15),rep(sqrt(.6),5))
P2 <- cbind(p1g2,p2g2,p3g2,p4g2)
PSI1 <- diag(1-rowSums(P1^2))
SIGMA1 <- P1%*%t(P1)+PSI1 
PSI2 <- diag(1-rowSums(P2^2))
SIGMA2 <- P2%*%t(P2)+PSI2
for (i in 1:length(n)){
  getValue <- function(){
    data1 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA1, empirical = TRUE)
    data2 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA2, empirical = TRUE)
    data3 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA1, empirical = TRUE)
    data4 <- mvrnorm(n = n[i], mu = rep(0,20), Sigma = SIGMA2, empirical = TRUE)
    DATA <- list(data1,data2,data3,data4)
    
    G <- 4
    R <- 4
    
    #model selection
    selmodel <- seqstrategy(DATA,R, SP=0.9, ALPH=0.1)
    
    #calculation performance measures
    parameters <- performancemeasures(selmodel,P1,P2,R,G)[[1]] 
    #parameters
    #round(cbind(P1,selmodel$loadings[[1]],P2,selmodel$loadings[[2]]),2)
    return(parameters)
  }
  
  result38 <- as.data.frame(t(replicate(50,getValue())))
  result38$size <- c(rep(n[i],50))
  result38$nrpcs <- 4
  result38$nrgroups <- 4
  final_result38 <- rbind(final_result38,result38)
}

#colnames(final_result38) <- c("Recovery","Tucker","Group_size")
save(final_result38, file = "final_result38v2.Rda")
