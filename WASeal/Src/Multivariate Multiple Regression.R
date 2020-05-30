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


dataCLIM <- subset(data,Year>=1960&TPear<=2008 & AA=="PRO")
dataCLIM <-data3  %>% select(MEI, 
                               PDO,
                               NPGO,
                               WA.SST.Su, 
                               UpInAn.45.Summer,
                               Col.Dis.high,
                               UpInAn.45.Spring,
                               TP.norm,
                               TP,
                              Year,
                               AA,
                               Sample.ID,
                               PHE.mean,
                               d13C.s,
                               Location.2)

dataCLIM <- dataCLIM[complete.cases(dataCLIM), ]

  aic.output <- rbind(extractAIC(lm(TP~Location.2, data=dataCLIM))[2], 
                      extractAIC(lm(TP~PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(TP~PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(TP~PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(TP~UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(TP~MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(TP~MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      extractAIC(lm(TP~MEI+NPGO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#7
                      extractAIC(lm(TP~PDO+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#8
                      extractAIC(lm(TP~PDO+MEI+NPGO+Location.2, data=dataCLIM))[2],#9
                      extractAIC(lm(TP~PDO+NPGO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      
                      extractAIC(lm(TP~Location.2+WA.SST.Su, data=dataCLIM))[2], 
                      extractAIC(lm(TP~WA.SST.Su+PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~WA.SST.Su+NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~WA.SST.Su+MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~WA.SST.Su+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~WA.SST.Su+PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(TP~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(TP~WA.SST.Su+PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(TP~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(TP~WA.SST.Su+MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(TP~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      
                      extractAIC(lm(TP~Location.2+UpInAn.45.Summer, data=dataCLIM))[2], 
                      extractAIC(lm(TP~UpInAn.45.Summer+PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~UpInAn.45.Summer+NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~UpInAn.45.Summer+MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~UpInAn.45.Summer+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~UpInAn.45.Summer+PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(TP~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(TP~UpInAn.45.Summer+PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(TP~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(TP~UpInAn.45.Summer+MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(TP~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      
                      extractAIC(lm(TP~Location.2+Col.Dis.high, data=dataCLIM))[2], 
                      extractAIC(lm(TP~Col.Dis.high+PDO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~Col.Dis.high+NPGO+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~Col.Dis.high+MEI+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~Col.Dis.high+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],
                      extractAIC(lm(TP~Col.Dis.high+PDO+NPGO+Location.2, data=dataCLIM))[2],#1
                      extractAIC(lm(TP~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#2
                      extractAIC(lm(TP~Col.Dis.high+PDO+MEI+Location.2, data=dataCLIM))[2],#3
                      extractAIC(lm(TP~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(TP~Col.Dis.high+MEI+NPGO+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(TP~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2, data=dataCLIM))[2],#6
                      
                      extractAIC(lm(TP~Col.Dis.high+WA.SST.Su+Location.2, data=dataCLIM))[2],#4
                      extractAIC(lm(TP~Col.Dis.high+UpInAn.45.Summer+Location.2, data=dataCLIM))[2],#5
                      extractAIC(lm(TP~WA.SST.Su+UpInAn.45.Summer+Location.2, data=dataCLIM))[2],#6
                      extractAIC(lm(TP~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2, data=dataCLIM))[2]#6
                      
                      
                      
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

########################     Prey Models       ############################



data5<- data3 %>% select(allSmolt, 
                         HakeBiomass,
                         Herring.Biomass,
                         Chinook, 
                         HarborSeal,
                         TP.ALA,
                        # TP.ASP,
                         TP.PRO,
                        #TP.VAL,
                         TP, #TP for GLU
                        Year,
                         Sample.ID, 
                         Location.2,
                       Chum,
                       Coho)

Prey <- data.frame(data5[complete.cases(data5), ])



Prey.aic.output <- rbind(extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+Chinook+HakeBiomass+Location.2, data=Prey))[2],#7
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#8
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+allSmolt+Chinook+Location.2, data=Prey))[2],#9
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+Chinook+HakeBiomass+Location.2, data=Prey))[2],
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2+Chum, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                     # extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2+Coho, data=Prey))[2], 
                    #  extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+Location.2, data=Prey))[2],
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Chinook+Location.2, data=Prey))[2],
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+allSmolt+Location.2, data=Prey))[2],
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+HakeBiomass+Location.2, data=Prey))[2],
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      #extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2+HarborSeal, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Chum+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Coho+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Coho+Location.2, data=Prey))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Coho+HarborSeal+Location.2, data=Prey))[2]#6
                      
                      
                      
  )
  
  model.names <- c("Location",
                   "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
                   "5. Chinook, Hatch", "6. Hatch, Hake", "7. Chinook, Hake, Hatch", "8. Herring, Hatch, Hake", "9. Herring, Hatch, Chinook", 
                   "11. Herring Hake Chinook",
                   
                   "Chum, Location", "Chum, Herring", "Chum, Chinook", "Chum, Hatch", "Chum, Hake","1. Chum,  Herring, Chinook", "2. Chum, Herring, Hake", 
                   "3. Chum, Herring, Hatch", "4. Chum, Chinook, Hake","5. Chum, Chinook, Hatch", "6. Chum, Hatch, Hake",
                   
                  # "Coho, Location", "Coho, Herring", "Coho, Chinook", "Coho, Hatch", "Coho, Hake","1. Coho, Herring, Chinook", "2.Coho,  Herring, Hake",
                 #  "3. Coho, Herring, Hatch", "4.Coho,  Chinook, Hake","5.Coho,  Chinook, Hatch", "6. Coho, Hatch, Hake",
                   
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
fit6 <- lmer(mvbind(TP, TP.ALA, TP.PRO) ~ Herring.Biomass+(1|Location.2), data = Prey)


summarTP(fit6)
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



