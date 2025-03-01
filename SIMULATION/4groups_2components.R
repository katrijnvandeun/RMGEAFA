#4 groups and 2 components
#primary loading shifts
#no of loading differences = 4
#generated loading matrices are assigned to 2 groups
#setwd("C:/Users/ASUS/surfdrive/Traineeship 1 - Tra Le/Final result dataframes")
setwd("C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/cSPCA/R/Simulation")
library(MASS)
n <- c(50,200,600,1000)
final_result21 <- data.frame()
set.seed(8790)
p1g1 <- c(0,rep(sqrt(.6),9),rep(0,10))
p2g1 <- c(sqrt(.6),rep(0,9),rep(sqrt(.6),10))
P1 <- cbind(p1g1,p2g1)
p1g2 <- c(rep(sqrt(.6),2),0,rep(sqrt(.6),7),rep(0,10))
p2g2 <- c(0,0,sqrt(.6),rep(0,7),rep(sqrt(.6),10))
P2 <- cbind(p1g2,p2g2)
#p1g3 <- c(rep(sqrt(.6),10),rep(0,10))
#p2g3 <- c(rep(0,10),rep(.6),10)
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
  result21 <- as.data.frame(t(replicate(50,getValue())))
  result21$size <- c(rep(n[i],50))
  result21$nrpcs <- 2
  result21$nrgroups <- 4
  final_result21 <- rbind(final_result21,result21)
}

save(final_result21, file = "final_result21v2.Rda")
