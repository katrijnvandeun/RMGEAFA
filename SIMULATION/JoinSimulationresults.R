######## Simulation Study JSPCA #################
### Codes to generate boxplots
## Author: Tra Le & Katrijn Van Deun

#####################################################################################

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
RESULT <- c()
for (i in 1:40){
  load(paste0("final_result",i,"v2.Rda"))
}

final_result1 <- final_result
# one big dataframe for 2 groups 2 components
dataG2R2 <- rbind(final_result1, final_result2, final_result3, final_result4,
              final_result5, final_result6, final_result7, final_result8,
              final_result9, final_result10)
dataG2R2$condition <- c(rep("PL shift", 400), rep("CL0.4", 200), rep("CL0.2",200),
                    rep("CL0.4", 200), rep("CL0.2",200), rep("PL0.4", 200), rep("PL0.2", 200),
                    rep("PL0.4", 200), rep("PL0.2", 200))
dataG2R2$diff <- c(rep(4, 200), rep(16, 200), rep(4, 400), rep(16, 400), rep(4, 400), rep(16, 400))
RESULT <- rbind(RESULT,dataG2R2)

# one big dataframe for 2 groups 4 components
dataG2R4 <- rbind(final_result11, final_result12, final_result13, final_result14,
                  final_result15, final_result16, final_result17, final_result18,
                  final_result19, final_result20)
dataG2R4$condition <- c(rep("PL shift", 400), rep("CL0.4", 200), rep("CL0.2",200),
                        rep("CL0.4", 200), rep("CL0.2",200), rep("PL0.4", 200), rep("PL0.2", 200),
                        rep("PL0.4", 200), rep("PL0.2", 200))
dataG2R4$diff <- c(rep(4, 200), rep(16, 200), rep(4, 400), rep(16, 400), rep(4, 400), rep(16, 400))
RESULT <- rbind(RESULT,dataG2R4)

# one big dataframe for 4 groups 2 components
dataG4R2 <- rbind(final_result21, final_result22, final_result23, final_result24,
                  final_result25, final_result26, final_result27, final_result28,
                  final_result29, final_result30)
dataG4R2$condition <- c(rep("PL shift", 400), rep("CL0.4", 200), rep("CL0.2",200),
                        rep("CL0.4", 200), rep("CL0.2",200), rep("PL0.4", 200), rep("PL0.2", 200),
                        rep("PL0.4", 200), rep("PL0.2", 200))
dataG4R2$diff <- c(rep(4, 200), rep(16, 200), rep(4, 400), rep(16, 400), rep(4, 400), rep(16, 400))
RESULT <- rbind(RESULT,dataG4R2)

# one big dataframe for 4 groups 4 components 
colnames(final_result32)<-colnames(final_result13)#result32 mislabeled
dataG4R4 <- rbind(final_result31, final_result32, final_result33, final_result34,
                  final_result35, final_result36, final_result37, final_result38,
                  final_result39, final_result40)
dataG4R4$condition <- c(rep("PL shift", 400), rep("CL0.4", 200), rep("CL0.2",200),
                        rep("CL0.4", 200), rep("CL0.2",200), rep("PL0.4", 200), rep("PL0.2", 200),
                        rep("PL0.4", 200), rep("PL0.2", 200))
dataG4R4$diff <- c(rep(4, 200), rep(16, 200), rep(4, 400), rep(16, 400), rep(4, 400), rep(16, 400))
RESULT <- rbind(RESULT,dataG4R4)

labels <- c('FP_PC1','FP_PC2','FN_PC1','FN_PC2',
            'FUSION_FP','FUSION_FN','Tucker',
            'n','R','G','Condition','Diff')
colnames(RESULT) <- labels
RESULT$G[7601:7800]<-4#mislabeled
colSums(RESULT[,1:6]==0)#preliminary exploration
par(mar=c(2,2,2,2))
hist(rowSums(RESULT[,1:6]))#find difficult datasets/conditions 
RESULT[rowSums(RESULT[,1:6])>20,]#16 loading differences, CL0.4/PL shift, 4G4R
sum(RESULT[,7]<0.7)#1 TUcker congruence
sum(RESULT[,7]<0.85)#22
sum(RESULT[,7]<0.95)#39

#make table to compare with De Roover's results
TABLE <- matrix(nrow=16,ncol=3)  #%datasets no FP/FN for fusion and simple structure, everything correct
colnames(TABLE) <- c('FUSION', 'SIMPLE S.', 'BOTH')

#find number of datasets 100%correct: OVERALL
fusion_ISF <- sum(rowSums(RESULT[,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[,1:6]==0)==6)#ISF
TABLE[16,] <- c(fusion_ISF,simple_ISF,both_ISF)/8000

#find number of datasets 100%correct: Per nr of components
Rvec <- c(2,4)
Gvec <- c(2,4)
for (r in 1:2){
  sel <- which(RESULT$R==Rvec[r])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[r,] <- c(fusion_ISF,simple_ISF,both_ISF)/4000
}
#find number of datasets 100%correct: Per sample size
Nvec <- c(50,200,600,1000)
for (n in 1:4){
  sel <- which(RESULT$n==Nvec[n])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[2+n,] <- c(fusion_ISF,simple_ISF,both_ISF)/2000
}
#find number of datasets 100%correct: Per group size
for (q in 1:2){
  sel <- which(RESULT$G==Gvec[q])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[6+q,] <- c(fusion_ISF,simple_ISF,both_ISF)/4000
}
#specific conditions
#1.Primary loading shifts: 4800 data sets
sel <- which(RESULT$Condition=='PL shift' | RESULT$Condition=='PL0.2' | RESULT$Condition=='PL0.4')
fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
TABLE[9,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
#2.remaining conditions
conds <- c('CL0.4','CL0.2','PL0.4','PL0.2')
for (c in 1:4){
  sel <- which(RESULT$Condition==conds[c])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[9+c,] <- c(fusion_ISF,simple_ISF,both_ISF)/length(sel)
}
#find number of datasets 100%correct: Per group size
Dvec <- c(4,16)
for (d in 1:2){
  sel <- which(RESULT$Diff==Dvec[d])
  fusion_ISF <- sum(rowSums(RESULT[sel,5:6]==0)==2)#ISF
  simple_ISF <- sum(rowSums(RESULT[sel,1:4]==0)==4)#ISF
  both_ISF <- sum(rowSums(RESULT[sel,1:6]==0)==6)#ISF
  TABLE[13+d,] <- c(fusion_ISF,simple_ISF,both_ISF)/4000
}
write.table(TABLE, file='simresultsNEWv2.txt',sep = '\t',dec = ",")





# set grid labels for the graph
library(ggplot2)
library(ggsci)
Nlabels <- c('50' = 'N = 50', '200' = "N = 200", '600' = "N = 600", '1000' = "N = 1000")
glabels <- c('4' = "4 differences", '16' = "16 differences")

#Tucker congruence
# boxplot for 2 groups 2 components
pdf("Tucker.pdf", width = 6, height = 5)
ggplot(data = dataG2R2[dataG2R2$size==50,], aes(x=condition, y = V14))+
  
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Tucker Congruence") + 
  ylim(.6, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12))
dev.off()


## Recovery Rate --------------------------------------------------------------------------

# boxplot for 2 groups 2 components
pdf("./Recovery Rate G2R2.pdf", width = 8, height = 10)
ggplot(data = dataG2R2, aes(x=condition, y = Recovery))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Recovery Rate") + 
  ylim(.6, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Group_size ~ diff, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                       diff = as_labeller(glabels)))
dev.off()

# boxplot for 2 groups 4 components
png("Recovery Rate G2R4.png", width = 3100, height = 2100, res = 300)
ggplot(data = dataG2R4, aes(x=condition, y = Recovery))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Recovery Rate") + 
  ylim(.6, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=12), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(diff~Group_size, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                    diff = as_labeller(glabels)))
dev.off()

# boxplot for 4 groups 2 components
png("Recovery Rate G4R2.png", width = 2100, height = 3100, res = 300)
ggplot(data = dataG4R2, aes(x=condition, y = Recovery))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Recovery Rate") + 
  ylim(.4, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Group_size ~ diff, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                    diff = as_labeller(glabels)))
dev.off()

# boxplot for 4 groups 4 components
png("Recovery Rate G4R4.png", width = 2100, height = 3100, res = 300)
ggplot(data = dataG4R4, aes(x=condition, y = Recovery))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Recovery Rate") + 
  ylim(.6, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Group_size ~ diff, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                    diff = as_labeller(glabels)))
dev.off()

## Tucker's Congruence --------------------------------------------------------------------------

# boxplot for 2 groups 2 components
png("Tucker congruence G4R2.png", width = 2100, height = 3100, res = 300)
ggplot(data = dataG2R2, aes(x=condition, y = Tucker))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Tucker's Congruence") + 
  ylim(.6, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Group_size ~ diff, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                    diff = as_labeller(glabels)))
dev.off()

# boxplot for 2 groups 4 components
png("Tucker congruence G2R4.png", width = 2100, height = 3100, res = 300)
ggplot(data = dataG2R4, aes(x=condition, y = Tucker))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Tucker's Congruence") + 
  ylim(.6, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=12), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Group_size ~ diff, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                  diff = as_labeller(glabels)))
dev.off()

# boxplot for 4 groups 2 components
png("Tucker congruence G4R2.png", width = 2100, height = 3100, res = 300)
ggplot(data = dataG4R2, aes(x=condition, y = Tucker))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Tucker's Congruence") +  
  ylim(.4, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Group_size ~ diff, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                    diff = as_labeller(glabels)))
dev.off()

# boxplot for 4 groups 4 components
png("Tucker congruence G4R4.png", width = 2100, height = 3100, res = 300)
ggplot(data = dataG4R4, aes(x=condition, y = Tucker))+
  geom_boxplot()+
  xlab("Types of Loading Differences") +
  ylab("Tucker's Congruence") +  
  ylim(.6, 1) +
  theme(panel.background = element_rect(fill = "grey96"), axis.text=element_text(size=12), 
        axis.title=element_text(size=15), legend.text = element_text(size=15), 
        legend.title = element_text(size=15), strip.text = element_text(size = 12)) +
  facet_grid(Group_size ~ diff, labeller = labeller(Group_size = as_labeller(Nlabels),
                                                    diff = as_labeller(glabels)))
dev.off()