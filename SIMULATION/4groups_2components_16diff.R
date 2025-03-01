#4 groups and 2 components
#primary loading shifts
#no of loading differences = 16
#setwd("C:/Users/ASUS/surfdrive/Traineeship 1 - Tra Le/Final result dataframes")
setwd("C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/cSPCA/R/Simulation")
library(MASS)
n <- c(50,200,600,1000)
final_result22 <- data.frame()
set.seed(45468)
p1g1 <- c(0,0,rep(sqrt(.6),10),rep(0,8))
p2g1 <- c(sqrt(.6),sqrt(.6),rep(0,10),rep(sqrt(.6),8))
P1 <- cbind(p1g1,p2g1)
p1g2 <- c(rep(sqrt(.6),2),0,0,rep(sqrt(.6),6),0,0,rep(sqrt(.6),2),rep(0,6))
p2g2 <- c(0,0,rep(sqrt(.6),2),rep(0,6),rep(sqrt(.6),2),0,0,rep(sqrt(.6),6))
P2 <- cbind(p1g2,p2g2)
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
    R <- 2
    
    #model selection
    selmodel <- seqstrategy(DATA,R, SP=0.9, ALPH=0.1)
    
    #calculation performance measures
    parameters <- performancemeasures(selmodel,P1,P2,R,G)[[1]]    
    #round(cbind(P1,selmodel$loadings[[1]],P2,selmodel$loadings[[2]]),2)
    return(parameters)
  }
  result22 <- as.data.frame(t(replicate(50,getValue())))
  result22$size <- c(rep(n[i],50))
  result22$nrpcs <- 2
  result22$nrgroups <- 4
  final_result22 <- rbind(final_result22,result22)
}

#colnames(final_result22) <- c("Recovery","Tucker","Group_size")
save(final_result22, file = "final_result22v2.Rda")
