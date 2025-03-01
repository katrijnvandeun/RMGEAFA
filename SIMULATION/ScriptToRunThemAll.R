rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("C:/Users/kvandeun/surfdrive/Manuscripts/InProgress/cSPCA/R/Simulation")

source("JointSPCA_finalizing.R")
source("IS_JSPCA.R")
source("seqstrategy.R")
source("performancemeasures.R")
source("InitialLoadings.R")#2nd version: Rotation towards simple structure
source("pstr.R")
source("2groups_2components.R")          #finalresult1
source("2groups_2components_16differences.R")          #finalresult2
source("2groups_2components_CL.R")          #finalresult3
source("2groups_2components_CL_2.R")      #finalresult4 
source("2groups_2components_CL_16differences.R")      #finalresult5 
source("2groups_2components_CL_16differences_2.R")   #finalresult6
source("2groups_2components_PLdecreases.R")#??  #finalresult7
source("2groups_2components_PLdecreases_2.R") #finalresult8
source("2groups_2components_PLdecreases_16diff.R") #finalresult9
source("2groups_2components_PLdecreases_16diff_2.R") #finalresult10

source("2groups_4components.R")          #finalresult11
source("2groups_4components_16diff.R")          #finalresult12
source("2groups_4components_CL.R")          #finalresult13
source("2groups_4components_CL_2.R")      #finalresult14 
source("2groups_4components_CL_16diff.R")      #finalresult15 
source("2groups_4components_CL_16diff_2.R")   #finalresult16
source("2groups_4components_PLdecreases.R")#??  #finalresult17
source("2groups_4components_PLdecreases_2.R") #finalresult18
source("2groups_4components_PLdecreases_16diff.R") #finalresult19
source("2groups_4components_PLdecreases_16diff_2.R") #finalresult20

###all done
source("4groups_2components.R")      #finalresult21
source("4groups_2components_16diff.R")  #finalresult22   
source("4groups_2components_CL.R")     #finalresult23    
source("4groups_2components_CL_2.R")    #finalresult24   
source("4groups_2components_CL_16diff.R")  #finalresult25
source("4groups_2components_CL_16diff_2.R") #finalresult26
source("4groups_2components_PL.R")     
source("4groups_2components_PL_16diff.R")  
source("4groups_2components_PL_16diff_2.R")
source("4groups_2components_PL_2.R")    

source("4groups_4components.R")     #finalresult31
source("4groups_4components_16diff.R")  #finalresult32   
source("4groups_4components_CL.R")     #finalresult33    
source("4groups_4components_CL_2.R")    #finalresult34   
source("4groups_4components_CL_16diff.R")  #finalresult35
source("4groups_4components_CL_16diff_2.R") #finalresult36
source("4groups_4components_PL.R") #finalresult37
source("4groups_4components_PL_2.R") #finalresult38
source("4groups_4components_PL_16diff.R")  #finalresult39
source("4groups_4components_PL_16diff_2.R")#finalresult40
   
