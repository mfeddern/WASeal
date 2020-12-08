rm(list = ls())
data <- read.csv("Data/Compiled/HierarchicalData2.csv")

#install.packages("brms")
library(dplyr)
library(brms)
library(tidyr)
library(dplyr)
library(car)
library(mgcViz)
library(AICcmodavg)
library(vapniks)
require(nlme)
require(mgcv)
library(qpcR)
library(dotwhisker)
library(ggpubr)
library(mvtnorm)
library(lme4)
library(corrplot)
library(MARSS)

data2<- subset(data, AA=='Glu')
#df <- tibble(x = 1:3, TP = 3:1)
data3 <- cbind(data2, TP.ALA.n = subset(data, AA=='ALA')$TP.norm, 
               TP.VAL.n = subset(data, AA=='VAL')$TP.norm, 
               TP.ASP.n = subset(data, AA=='ASP')$TP.norm, 
               TP.PRO.n = subset(data, AA=='PRO')$TP.norm, 
               TP.ALA = subset(data, AA=='ALA')$TP, 
               TP.VAL = subset(data, AA=='VAL')$TP, 
               TP.ASP = subset(data, AA=='ASP')$TP, 
               TP.PRO = subset(data, AA=='PRO')$TP)
########################     Climate Models       ############################


dataCLIM <-data3  %>% select(MEI, 
                               PDO,
                               NPGO,
                               WA.SST.Su, 
                               UpInAn.45.Summer,
                               Col.Dis.high,
                               UpInAn.45.Spring,
                               TP.norm,
                               TP,
                             TP.ALA,
                             TP.ASP,
                             TP.PRO,
                             TP.VAL,
                              Year,
                               AA,
                               Sample.ID,
                               PHE.mean,
                               d13C.s,
                               Location.2)

dataCLIM <- dataCLIM[complete.cases(dataCLIM), ]

  aic.output <- rbind(extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Location.2, data=dataCLIM))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~MEI+NPGO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#7
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PDO+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#8
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PDO+MEI+NPGO+Location.2, data=dataCLIM))[2],#9
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PDO+NPGO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Location.2+WA.SST.Su, data=dataCLIM))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Location.2+UpInAn.45.Summer, data=dataCLIM))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Location.2+Col.Dis.high, data=dataCLIM))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+WA.SST.Su+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Col.Dis.high+UpInAn.45.Summer+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+UpInAn.45.Summer+Location.2, data=dataCLIM))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2, data=dataCLIM))[2]#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "PDO", "NPGO", "MEI", "Upwelling (Sp)","1. PDO, NPGO", "2. PDO, Upwelling (Sp)", "3. PDO, MEI", "4. NPGO, Upwelling (Sp)",
                   "5. NPGO, MEI", "6. MEI, Upwelling (Sp)", "7. NPGO, Upwelling (Sp), MEI", "8. PDO, MEI, Upwelling (Sp)", "9. PDO, MEI, NPGO", 
                   "11. PDO Upwelling (Sp) NPGO",
                   
                   "WA.SST.Su, Location", "WA.SST.Su, PDO", "WA.SST.Su, NPGO", "WA.SST.Su, MEI", "WA.SST.Su, Upwelling (Sp)","1. WA.SST.Su,  PDO, NPGO", "2. WA.SST.Su, PDO, Upwelling (Sp)", 
                   "3. WA.SST.Su, PDO, MEI", "4. WA.SST.Su, NPGO, Upwelling (Sp)","5. WA.SST.Su, NPGO, MEI", "6. WA.SST.Su, MEI, Upwelling (Sp)",
                   
                   "UpInAn.45.Summer, Location", "UpInAn.45.Summer, PDO", "UpInAn.45.Summer, NPGO", "UpInAn.45.Summer, MEI", "UpInAn.45.Summer, Upwelling (Sp)","1. UpInAn.45.Summer, PDO, NPGO", "2.UpInAn.45.Summer,  PDO, Upwelling (Sp)",
                   "3. UpInAn.45.Summer, PDO, MEI", "4.UpInAn.45.Summer,  NPGO, Upwelling (Sp)","5.UpInAn.45.Summer,  NPGO, MEI", "6. UpInAn.45.Summer, MEI, Upwelling (Sp)",
                   
                   "Col.Dis.high, Location", "Col.Dis.high, PDO", "Col.Dis.high, NPGO", "Col.Dis.high, MEI", "Upwelling (Sp)","1. Col.Dis.high, PDO, NPGO",
                   "2.Col.Dis.high,  PDO, Upwelling (Sp)", "3.Col.Dis.high,  PDO, MEI", "4.Col.Dis.high,  NPGO, Upwelling (Sp)","5. Col.Dis.high, NPGO, MEI", 
                   "6. Col.Dis.high, MEI, Upwelling (Sp)",
                   
                   "Col.Dis.high, WA.SST.Su", " Col.Dis.high UpInAn.45.Summer", "WA.SST.Su, UpInAn.45.Summer", "WA.SST.Su, UpInAn.45.Summer, Col.Dis.high")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summarTP[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAIC", "AICc Weight")
  x<-data.frame(aic.output)
  subset(x, delAIC<=2)

  Clim.mod.1 <- lm(cbind(TP, TP.ALA, TP.PRO, TP.VAL) ~ PDO+MEI+WA.SST.Su+Location.2, data = dataCLIM)
  summary(Clim.mod.1)
acf(Clim.mod.1)
durbinWatsonTest(Clim.mod.1)
  ########################     Prey Models       ############################



data5<- data3 %>% select(allSmolt, 
                         HakeBiomass,
                         Herring.Biomass,
                         Chinook, 
                         HarborSeal,
                         TP.ALA,
                        TP.ASP,
                         TP.PRO,
                        TP.VAL,
                         TP, #TP for GLU
                        Year,
                         Sample.ID, 
                         Location.2,
                       Chum,
                       Coho)

Prey <- data.frame(data5[complete.cases(data5), ])



Prey.aic.output <- rbind(extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Location.2, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Herring.Biomass*HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~allSmolt+Chinook+HakeBiomass+Location.2, data=Prey))[2],#7
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Herring.Biomass*HakeBiomass+allSmolt+Location.2, data=Prey))[2],#8
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Herring.Biomass+allSmolt+Chinook+Location.2, data=Prey))[2],#9
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Herring.Biomass*HakeBiomass+Chinook+Location.2, data=Prey))[2],
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Location.2+Chum, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+Herring.Biomass*HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                     extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Location.2+Coho, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+Herring.Biomass*HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Coho+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Location.2+HarborSeal, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+Herring.Biomass*HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+Chum+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~HarborSeal+Coho+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+Coho+Location.2, data=Prey))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL) ~Chum+Coho+HarborSeal+Location.2, data=Prey))[2]#6
                      
                      
                      
  )
  
  model.names <- c("Location",
                   "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
                   "5. Chinook, Hatch", "6. Hatch, Hake", "7. Chinook, Hake, Hatch", "8. Herring, Hatch, Hake", "9. Herring, Hatch, Chinook", 
                   "11. Herring Hake Chinook",
                   
                   "Chum, Location", "Chum, Herring", "Chum, Chinook", "Chum, Hatch", "Chum, Hake","1. Chum,  Herring, Chinook", "2. Chum, Herring, Hake", 
                   "3. Chum, Herring, Hatch", "4. Chum, Chinook, Hake","5. Chum, Chinook, Hatch", "6. Chum, Hatch, Hake",
                   
                   "Coho, Location", "Coho, Herring", "Coho, Chinook", "Coho, Hatch", "Coho, Hake","1. Coho, Herring, Chinook", "2.Coho,  Herring, Hake",
                   "3. Coho, Herring, Hatch", "4.Coho,  Chinook, Hake","5.Coho,  Chinook, Hatch", "6. Coho, Hatch, Hake",
                   
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook", "Harbor Seal, Hatch", "Hake","1. Harbor Seal, Herring, Chinook",
                   "2.Harbor Seal,  Herring, Hake", "3.Harbor Seal,  Herring, Hatch", "4.Harbor Seal,  Chinook, Hake","5. Harbor Seal, Chinook, Hatch", 
                   "6. Harbor Seal, Hatch, Hake",
                   
                   "Harbor Seal, Chum", " Harbor Seal Coho", "Chum, Coho", "Chum, Coho, Harbor Seal")
  
  row.names(Prey.aic.output) <- model.names
  delaic <- Prey.aic.output-min(Prey.aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summarTP[,14]
  Prey.aic.output <- cbind(Prey.aic.output, delaic, aic.weight)
  
  colnames(Prey.aic.output)<- c("AIC", "delAIC", "AIC Weight")
  x<-data.frame(Prey.aic.output)

  subset(x, delAIC<=2)
Prey.mod.1 <- lm(cbind(TP, TP.ALA, TP.PRO,TP.VAL) ~ Herring.Biomass*HakeBiomass+allSmolt+Location.2, data = Prey)
acf(Prey.mod.1)

summary(fit6)
vcov(fit6)
Manova(fit6)
Anova(fit6)$aic
vif(fit6)



fit6 <- lm(cbind(TP, TP.ALA, TP.PRO) ~ Location.2, data = Prey)
summarTP(fit6)
vcov(fit6)
summarTP(Manova(fit6))
Anova(fit6)
extractAIC(fit6)




########################     Nutrient Models       ############################

#dataNutrient <- subset(data,  AA=="Glu")
dataNut <-data3%>% select(PHE.mean,
                                  PHE.norm,
                                  AA,
                                  d13C.s, 
                                  d13C.norm,
                                  TP.norm,
                                  TP,
                          TP.ALA,
                          TP.VAL,
                          TP.PRO,
                                  Year,
                                  Sample.ID,
                                  Location.2)
dataNut <- dataNut[complete.cases(dataNut), ]


  
aic.output <- rbind(
    extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~Location.2, data=dataNut))[2],
    extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PHE.norm+d13C.norm+Location.2, data=dataNut))[2],
    extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~PHE.norm+Location.2, data=dataNut))[2],
    extractAIC(lm(cbind(TP, TP.ALA, TP.PRO,   TP.VAL)~d13C.norm+Location.2, data=dataNut))[2])
  
  names <- seq(1,n,1)
  model.names <- c( "Location","Phe, 13C", "PHE", "13C")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")

  x<-data.frame(aic.output)
  
  subset(x, delAICc<=2)
  NUTR.mod.1<- lm(cbind(TP, TP.ALA, TP.PRO,TP.VAL) ~ PHE.norm+d13C.s++Location.2, data =dataNut )
summary(NUTR.mod.1)  





########################     Full Models       ############################


########## Predict! #########
Prey.TP <- lm(TP ~ HakeBiomass*Herring.Biomass+allSmolt+Location.2, data = Prey)
Prey.TP.ALA <- lm(TP.ALA ~ HakeBiomass*Herring.Biomass+allSmolt+Location.2, data = Prey)
Prey.TP.PRO <- lm(TP.PRO ~ HakeBiomass*Herring.Biomass+allSmolt+Location.2, data = Prey)
Prey.TP.VAL <- lm(TP.VAL ~ HakeBiomass*Herring.Biomass+allSmolt+Location.2, data = Prey)

new.DATA.low <- data.frame(
  TP=Prey$TP, TP.ALA=Prey$TP.ALA,TP.PRO=Prey$TP.PRO, TP.VAL=Prey$TP.VAL,
  HakeBiomass = Prey$HakeBiomass,
  Herring.Biomass = rep(-1.62, 93), allSmolt=Prey$allSmolt,
  Location.2=Prey$Location.2
)

new.DATA.high <- data.frame(
  TP=Prey$TP, TP.ALA=Prey$TP.ALA,TP.PRO=Prey$TP.PRO, TP.VAL=Prey$TP.VAL,
  HakeBiomass = Prey$HakeBiomass,
  Herring.Biomass = rep(1.5, 93), allSmolt=Prey$allSmolt,
  Location.2=Prey$Location.2
)

pred<- predict(Prey.mod.1, new.DATA.low, interval = "confidence")
pred.low <-data.frame(predict(Prey.mod.1, new.DATA.low, interval = "confidence"))
pred.low.mod<- data.frame(TP=c(pred.low[,1], pred.low[,2], pred.low[,3], pred.low[,4]))

pred.high <-data.frame(predict(Prey.mod.1, new.DATA.high, interval = "confidence"))
pred.high.mod<- data.frame(TP=c(pred.high[,1], pred.high[,2], pred.high[,3], pred.high[,4]))


PRED.low  <- pred.low.mod %>% 
  add_column("Model" = c(rep('Glu',93), rep('Ala',93), rep('Pro',93),rep('Val',93))) %>%  
  add_column("allSmolt"=rep(new.DATA.low$allSmolt, 4)) %>%  
  add_column("HakeBiomass"=rep(new.DATA.low$HakeBiomass, 4)) %>%  
  add_column("Location.2"=rep(new.DATA.low$Location.2, 4)) %>%  
  add_column("Type"=rep("Low", 372))%>%  
  add_column("lwr"=low[,2]) %>%  
  add_column("upr"=low[,3])   
 
PRED.high  <- pred.high.mod %>% 
  add_column("Model" = c(rep('Glu',93), rep('Ala',93), rep('Pro',93),rep('Val',93))) %>%  
  add_column("allSmolt"=rep(new.DATA.high$allSmolt, 4)) %>%  
  add_column("HakeBiomass"=rep(new.DATA.high$HakeBiomass, 4)) %>%  
  add_column("Location.2"=rep(new.DATA.high$Location.2, 4)) %>%  
  add_column("Type"=rep("High", 372)) 


color9<- rep(c('#CCA65A','#7EBA68','#6FB1E7','#D494E1'), each=93)
color8<- rep(c('#7EBA68','#CCA65A','#6FB1E7','#D494E1'), each=80)
legend.size<-12

Herring.low.plot<-ggplot(PRED.low, aes(x=HakeBiomass, y=TP, color = Model, group = Model)) + 
  geom_point() + 
  geom_smooth(method=lm, alpha = .25, aes(fill = Model), se=F)+
  scale_colour_manual(values=c('#7EBA68','#CCA65A','#6FB1E7','#D494E1'))+
  theme_bw() + xlab("Hake Biomass") + ylab("Trophic Position") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Herring", subtitle="Low Biomass")+
  scale_y_continuous(breaks=c(2,4,6,8), limits=c(2,8))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0.5, 0.5, 0.5, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = 0,
        legend.title=element_text(size=legend.size),
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        legend.position = "none",
        axis.title.x=element_text(size=13),
        axis.title.y=element_text(size=13),
        
        plot.subtitle = element_text(hjust = 0.5, size=13))

Herring.high.plot <- ggplot(PRED.high, aes(x=HakeBiomass, y=TP, color = Model, group = Model)) + 
  geom_point() + 
  geom_smooth(method=lm, alpha = .25, aes(fill = Model), se=F)+
  scale_colour_manual(values=c('#7EBA68','#CCA65A','#6FB1E7','#D494E1'))+
  theme_bw() + xlab("Hake Biomass") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Herring", subtitle="High Biomass")+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0.5, 0.5, 0.5, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = 0,
        legend.title=element_text(size=legend.size),
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        axis.title.x=element_text(size=13),
        axis.title.y=element_text(size=13),
        plot.subtitle = element_text(hjust = 0.5, size=13))

########## PLOTS! #########
library(dotwhisker)
library(broom)
library(dplyr)
library(colorspace)
library(ggpubr)

color3<-  rep(c('#CCA65A','#7EBA68','#6FB1E7','#D494E1'), each=5)
color4<- rep(c('#CCA65A','#7EBA68','#6FB1E7','#D494E1'), times=5)


ENV <- tidy(Clim.mod.1)
ENV  <- ENV  %>% 
  add_column("model" = c(rep('Glu',5), rep('Ala',5), rep('Pro',5),rep('Val',5)))%>% 
  rename(labs = term)%>% 
  add_column("term" = rep(c('Intercept', "PDO", "MEI", 'Summer SST', 'Salish Sea'),4))


ENV.plot<-small_multiple(ENV, model='response') +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("Coefficient Estimate") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Environment") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1))


color3<- rep(c('#CCA65A','#7EBA68','#6FB1E7','#D494E1'), each=6)
color4<- rep(c('#CCA65A','#7EBA68','#6FB1E7','#D494E1'), times=6)


PREY <- tidy(Prey.mod.1)
PREY<-PREY  %>% 
  add_column("model" = c(rep('Glu',6), rep('Ala',6), rep('Pro',6),rep('Val',6)))%>% 
  rename(labs = term)%>% 
  add_column("term" = rep(c('Intercept', "Herring", "Hake", 'Smolts', 'Salish Sea', "Herring:Hake"),4))

PREY.plot<- small_multiple(PREY, model='response') +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1))


color3<- rep(c('#CCA65A','#7EBA68','#6FB1E7','#D494E1'), each=4)
color4<- rep(c('#CCA65A','#7EBA68','#6FB1E7','#D494E1'), times=4)

NUTR <- tidy(NUTR.mod.1)
NUTR<-NUTR  %>% 
  add_column("model" = c(rep('Glu',4), rep('Ala',4), rep('Pro',4),rep('Val',4)))%>% 
  rename(labs = term)%>% 
  add_column("term" = rep(c('Intercept', "Phenylalanine", "d13C", 'Salish Sea'),4))

NUTR.plot<-small_multiple(NUTR, model='response') +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Primary Productivity") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1))


sjPlot::plot_model(Prey.mod.1, type = "int", col=c(col[1], col[2]))+
  theme_bw() + xlab("Hake Biomass") + ylab("Trophic Position") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web", subtitle="Interactions")+
  labs(fill="Herring Biomass")+
  # scale_color_manual(values= c("#CCA65A", "#7EBA68"),name="Herring Biomass", labels = c("Low", "High"))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0.5, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = 0,
        legend.title=element_text(size=legend.size),
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))






pdf(file="Results/Presentation Figures/NutrientDRAFT.pdf", width=5, height=5)
NUTR.plot
dev.off()

pdf(file="Results/Presentation Figures/PreyDRAFT.pdf", width=5, height=5)
PREY.plot
dev.off()

pdf(file="Results/Presentation Figures/EnvDRAFT.pdf", width=5, height=5)
ENV.plot
dev.off()

pdf(file="Results/Presentation Figures/FullDRAFT.pdf", width=5, height=5)
FULL.plot
dev.off()

pdf(file="Results/Presentation Figures/HerringLow.pdf", width=6, height=5)
Herring.low.plot
dev.off()

pdf(file="Results/Presentation Figures/HerringHigh.pdf", width=6, height=5)
Herring.high.plot
dev.off()


pdf(file="Results/Figures/CoefPlot2.pdf", width=12, height=6)
ggarrange(ENV.plot, NUTR.plot, PREY.plot, rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 1, align= 'hv',
          heights =c(5))

dev.off()

pdf(file="Results/Figures/InteractPlot2.pdf", width=9, height=4)
ggarrange(Herring.low.plot, Herring.high.plot, rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 1, align= 'hv',
          heights =c(5), widths=c(5.25,6))

dev.off()
