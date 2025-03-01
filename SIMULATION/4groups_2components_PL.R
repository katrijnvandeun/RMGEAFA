#4 groups and 2 components
#primary loading decrease 0.4
#no of loading differences = 4
#setwd("C:/Users/ASUS/surfdrive/Traineeship 1 - Tra Le/Final result dataframes")
setwd("C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/cSPCA/R/Simulation")
library(MASS)
n <- c(50,200,600,1000)
final_result27 <- data.frame()
set.seed(8790)
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
  
  result27 <- as.data.frame(t(replicate(50,getValue())))
  result27$size <- c(rep(n[i],50))
  result27$nrpcs <- 2
  result27$nrgroups <- 4
  final_result27 <- rbind(final_result27,result27)
}

#colnames(final_result27) <- c("Recovery","Tucker","Group_size")
save(final_result27, file = "final_result27v2.Rda")
