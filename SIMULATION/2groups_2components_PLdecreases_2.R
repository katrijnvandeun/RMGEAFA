#2 groups and 2 components
#primary loading decrease 0.2
#no of loading differences = 4
library(MASS)
n <- c(50,200,600,1000)
final_result8 <- data.frame()
set.seed(45468)
p1g1 <- c(sqrt(.6)-.2,rep(sqrt(.6),9),rep(0,10))
p2g1 <- c(rep(0,10),sqrt(.6)-.2,rep(sqrt(.6),9))
P1 <- cbind(p1g1,p2g1)
p1g2 <- c(rep(sqrt(.6),4),sqrt(.6)-.2,rep(sqrt(.6),5),rep(0,10))
p2g2 <- c(rep(0,10),rep(sqrt(.6),4),sqrt(.6)-.2,rep(sqrt(.6),5))
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
    modelsel <- seqstrategy(DATA,R, SP=0.9, ALPH=0.1, TOL=3)
    
    #calculation performance measures
    parameters <- performancemeasures(modelsel$selmodel,P1,P2,R,G)[[1]]    
    #round(cbind(P1,modelsel$selmodel$loadings[[1]],P2,modelsel$selmodel$loadings[[2]]),2)
    return(parameters)
  }
  
  result8 <- as.data.frame(t(replicate(50,getValue())))
  result8$size <- c(rep(n[i],50))
  result8$nrpcs <- 2
  result8$nrgroups <- 2
  final_result8 <- rbind(final_result8,result8)
}

#colnames(final_result8) <- c("Recovery","Tucker","Group_size")
save(final_result8, file = "final_result8v2.Rda")