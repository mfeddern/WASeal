rm(list = ls())

library(devtools)
install_version("lme4", "1.1-25")
library(lme4)
library(AICcmodavg)

dataALL <-  read.csv("Data/Compiled/HierarchicalData.csv")
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag1<-data

dataALL <-  read.csv("Data/Compiled/HierarchicalData2.csv")
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag2<-data

dataALL <-  read.csv("Data/Compiled/HierarchicalData3.csv")
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag3<-data

dataALL <-  read.csv("Data/Compiled/HierarchicalData4.csv")
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag4<-data


dataALL <-  read.csv("Data/Compiled/HierarchicalData5.csv")
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag5<-data
########################    Hier Clim2 Models       ############################

dataCLIM <- subset(data.lag1, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataCLIM.1 <-dataCLIM %>% select(MEI, 
                               PDO,
                               NPGO,
                               WA.SST.Su, 
                               UpInAn.45.Spring,
                               Col.Dis.high,
                               UpInAn.45.Summer,
                               TP,
                               Year,
                               AA,
                               Sample.ID,
                               Location.2)
dataCLIM.1 <- dataCLIM.1[complete.cases(dataCLIM.1), ]

dataCLIM <- subset(data.lag2, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataCLIM.2 <-dataCLIM %>% select(MEI, 
                                   PDO,
                                   NPGO,
                                   WA.SST.Su, 
                                   UpInAn.45.Spring,
                                   Col.Dis.high,
                                   UpInAn.45.Summer,
                                   TP,
                                   Year,
                                   AA,
                                   Sample.ID,
                                   Location.2)
dataCLIM.2 <- dataCLIM.2[complete.cases(dataCLIM.2), ]

dataCLIM<- subset(data.lag3, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataCLIM.3 <-dataCLIM %>% select(MEI, 
                                   PDO,
                                   NPGO,
                                   WA.SST.Su, 
                                   UpInAn.45.Spring,
                                   Col.Dis.high,
                                   UpInAn.45.Summer,
                                   TP,
                                   Year,
                                   AA,
                                   Sample.ID,
                                   Location.2)
dataCLIM.3 <- dataCLIM.3[complete.cases(dataCLIM.3), ]

dataCLIM<- subset(data.lag4, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataCLIM.4 <-dataCLIM %>% select(MEI, 
                                 PDO,
                                 NPGO,
                                 WA.SST.Su, 
                                 UpInAn.45.Spring,
                                 Col.Dis.high,
                                 UpInAn.45.Summer,
                                 TP,
                                 Year,
                                 AA,
                                 Sample.ID,
                                 Location.2)
dataCLIM.4 <- dataCLIM.4[complete.cases(dataCLIM.4), ]

ModelSelection.CLIM<- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe)), #1
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe)), #2
                      AICc(lmer(y~PDO+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~MEI+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      AICc(lmer(y~PDO+NPGO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#7
                      #AICc(lmer(y~PDO+MEI+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#8
                      #AICc(lmer(y~MEI+NPGO+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#9
                      #AICc(lmer(y~MEI+NPGO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#7
                      #AICc(lmer(y~PDO+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#8
                      #AICc(lmer(y~PDO+MEI+NPGO+Location.2+(1|AA), data=dataframe)),#9
                      #AICc(lmer(y~PDO+NPGO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      
                      AICc(lmer(y~Location.2+(1|AA)+WA.SST.Su, data=dataframe)), #10
                      AICc(lmer(y~WA.SST.Su+PDO+Location.2+(1|AA), data=dataframe)),#11
                      AICc(lmer(y~WA.SST.Su+NPGO+Location.2+(1|AA), data=dataframe)),#12
                      AICc(lmer(y~WA.SST.Su+MEI+Location.2+(1|AA), data=dataframe)),#13
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#14
                      # AIC(lmer(y~WA.SST.Su+PDO+NPGO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#15
                      #AIC(lmer(y~WA.SST.Su+PDO+MEI+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#16
                      #AIC(lmer(y~WA.SST.Su+MEI+NPGO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#17
                      AICc(lmer(y~Location.2+(1|AA)+UpInAn.45.Summer, data=dataframe)), #18
                      AICc(lmer(y~UpInAn.45.Summer+PDO+Location.2+(1|AA), data=dataframe)),#19
                      AICc(lmer(y~UpInAn.45.Summer+NPGO+Location.2+(1|AA), data=dataframe)),#20
                      AICc(lmer(y~UpInAn.45.Summer+MEI+Location.2+(1|AA), data=dataframe)),#21
                      #AIC(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#22
                      # AIC(lmer(y~UpInAn.45.Summer+PDO+NPGO+Location.2+(1|AA), data=dataframe)),
                      #AICc(lmer(y~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      #AIC(lmer(y~UpInAn.45.Summer+PDO+MEI+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#22
                      #AIC(lmer(y~UpInAn.45.Summer+MEI+NPGO+Location.2+(1|AA), data=dataframe)),
                      #AICc(lmer(y~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#24
                      
                      AICc(lmer(y~Location.2+(1|AA)+Col.Dis.high, data=dataframe)), #23
                      AICc(lmer(y~Col.Dis.high+PDO+Location.2+(1|AA), data=dataframe)),#24
                      AICc(lmer(y~Col.Dis.high+NPGO+Location.2+(1|AA), data=dataframe)),#25
                      AICc(lmer(y~Col.Dis.high+MEI+Location.2+(1|AA), data=dataframe)),#26
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#27
                      # AICc(lmer(y~Col.Dis.high+PDO+NPGO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#28
                      #AICc(lmer(y~Col.Dis.high+PDO+MEI+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#29
                      # AICc(lmer(y~Col.Dis.high+MEI+NPGO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#30
                      
                      AICc(lmer(y~Col.Dis.high+WA.SST.Su+Location.2+(1|AA), data=dataframe)),#31
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#32
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#33
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2+(1|AA), data=dataframe))

                      
                      
                      
  )
  
  model.names <- c("1. Null","2. Location", "3. PDO", "4. NPGO", "5. MEI", "6. Upwelling (Sp)","NPGO,PDO", "7. PDO, Upwelling (Sp)",  
                   "8. NPGO, Upwelling (Sp)","9. MEI, Upwelling (Sp)",  
                   
                   "10. WA.SST.Su, Location", "11. WA.SST.Su, PDO", "12. WA.SST.Su, NPGO", "13. WA.SST.Su, MEI",
                   "14. WA.SST.Su, Upwelling (Sp)", "15. WA.SST.Su, PDO, Upwelling (Sp)", 
                   "16. WA.SST.Su, NPGO, Upwelling (Sp)", "17. WA.SST.Su, MEI, Upwelling (Sp)",
                   
                   "18. UpInAn.45.Summer, Location", "19. UpInAn.45.Summer, PDO", "20. UpInAn.45.Summer, NPGO", 
                   "21. UpInAn.45.Summer, MEI", #"22. UpInAn.45.Summer, Upwelling (Sp)",
                   #"23.UpInAn.45.Summer,  PDO, Upwelling (Sp)",
                   "22. UpInAn.45.Summer,  NPGO, Upwelling (Sp)", #"25. UpInAn.45.Summer, MEI, Upwelling (Sp)",
                   
                   "23. Col.Dis.high, Location", "24. Col.Dis.high, PDO", "25. Col.Dis.high, NPGO", "26. Col.Dis.high, MEI", "27. Upwelling (Sp)",
                   "28.Col.Dis.high,  PDO, Upwelling (Sp)", "29.Col.Dis.high,  NPGO, Upwelling (Sp)",
                   "30. Col.Dis.high, MEI, Upwelling (Sp)",
                   
                   "31. Col.Dis.high, WA.SST.Su", "32. Col.Dis.high UpInAn.45.Summer", 
                   "33. WA.SST.Su, UpInAn.45.Summer", "34. WA.SST.Su, UpInAn.45.Summer, Col.Dis.high")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}

model.selectionCLIM.1 <- ModelSelection.CLIM(dataCLIM.1, n, dataCLIM.1$TP)
clim.1 <-data.frame(model.selectionCLIM.1)
subset(clim.1, delAICc<=1.97)
modelENV.1 <- lmer(TP~Location.2+UpInAn.45.Summer+(1|AA), data=dataCLIM.1)
summary(modelENV.1)
modelENV.1 <- lmer(TP~Location.2+UpInAn.45.Summer+(1|AA), data=dataCLIM.1)
summary(lmer(TP~(1|AA), data=dataCLIM.1))

model.selectionCLIM.2 <- ModelSelection.CLIM(dataCLIM.2, n, dataCLIM.2$TP)
clim.2 <-data.frame(model.selectionCLIM.2)
subset(clim.2, delAICc<=1.97)
modelENV.2 <- lmer(TP~Location.2+WA.SST.Su+(1|AA), data=dataCLIM.2)
summary(modelENV.2)

model.selectionCLIM.3 <- ModelSelection.CLIM(dataCLIM.3, n, dataCLIM.3$TP)
clim.3 <-data.frame(model.selectionCLIM.3)
subset(clim.3, delAICc<=2)
modelENV.3 <- lmer(TP~Location.2+UpInAn.45.Summer+MEI+(1|AA), data=dataCLIM.3)
summary(modelENV.3)

model.selectionCLIM.4 <- ModelSelection.CLIM(dataCLIM.4, n, dataCLIM.4$TP)
clim.4 <-data.frame(model.selectionCLIM.4)
subset(clim.4, delAICc<=2)
modelENV.4 <- lmer(TP~Location.2+WA.SST.Su+(1|AA), data=dataCLIM.4)
summary(modelENV.4)

########################    Hier PREY Models       ############################

dataPrey <- subset(data.lag1, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataPrey.1 <-dataPrey %>% select(allSmolt, 
                            HakeBiomass,
                            Herring.Biomass,
                            AA,
                            Chinook, 
                            HarborSeal,
                            WildProduction,
                            HatcherySmolts,
                            Chum,
                            MEI,
                            TP,
                            Year,
                            Sample.ID,
                            Location.2,
                            eq,
                            beta)
dataPrey.1 <- dataPrey.1[complete.cases(dataPrey.1), ]

dataPrey <- subset(data.lag2, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataPrey.2 <-dataPrey %>% select(allSmolt, 
                                  HakeBiomass,
                                  Herring.Biomass,
                                  AA,
                                  Chinook, 
                                  HarborSeal,
                                 WildProduction,
                                 HatcherySmolts,
                                  MEI,
                                  TP,
                                  Year,
                                  Sample.ID,
                                  Location.2,
                                  eq,
                                  beta)
dataPrey.2 <- dataPrey.2[complete.cases(dataPrey.2), ]

dataPrey <- subset(data.lag3, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataPrey.3 <-dataPrey %>% select(allSmolt, 
                                  HakeBiomass,
                                  Herring.Biomass,
                                  AA,
                                  Chinook, 
                                  HarborSeal,
                                 WildProduction,
                                 HatcherySmolts,
                                  TP,
                                  Year,
                                  Sample.ID,
                                  Location.2,
                                  eq,
                                  beta)
dataPrey.3 <- dataPrey.3[complete.cases(dataPrey.3), ]

dataPrey <- subset(data.lag4, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataPrey.4 <-dataPrey %>% select(allSmolt, 
                                 HakeBiomass,
                                 Herring.Biomass,
                                 AA,
                                 Chinook, 
                                 HarborSeal,
                                 WildProduction,
                                 HatcherySmolts,
                                 TP,
                                 Year,
                                 Sample.ID,
                                 Location.2,
                                 eq,
                                 beta)
dataPrey.4 <- dataPrey.4[complete.cases(dataPrey.4), ]

dataPrey <- subset(data.lag4, AA=="GLU"|AA=="ALA"|AA=="VAL"|AA=="PRO")
dataPrey.5 <-dataPrey %>% select(allSmolt, 
                                 HakeBiomass,
                                 Herring.Biomass,
                                 AA,
                                 Chinook, 
                                 HarborSeal,
                                 WildProduction,
                                 HatcherySmolts,
                                 TP,
                                 Year,
                                 Sample.ID,
                                 Location.2,
                                 eq,
                                 beta)
dataPrey.5 <- dataPrey.5[complete.cases(dataPrey.5), ]

model.selection.PREY <- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe)), 
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Location.2+WildProduction+(1|AA), data=dataframe)),
                      AICc(lmer(y~Location.2+HatcherySmolts+(1|AA), data=dataframe)),
                      AICc(lmer(y~Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Chinook+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~allSmolt+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),#7
                      AICc(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#8
                     AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#9
                      AICc(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      
             
                      AICc(lmer(y~HarborSeal+Location.2+(1|AA), data=dataframe)), 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Chinook+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+allSmolt+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      #AICc(lmer(y~HarborSeal+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                
                  
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
                     AICc(lmer(y~Herring.Biomass*HakeBiomass+allSmolt+Location.2+(1|AA), data=dataframe))
                      
                      
                      
  )
  
  model.names <- c("Null","Location","Wild Smolts", "Hatch Smolts", "Herring", "Chinook", "Smolts", "Hake","1. Herring, Chinook", 
                   "2. Herring, Hake", "3. Herring, Smolts", "4. Chinook, Hake","CH SM",
                   "6. Smolts, Hake", "CH SM", "8. Herring, Smolts, Hake","CH SM",
                   "11. Herring Hake Chinook",
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook", "Harbor Seal, Smolts", 
                   "Hake","1. Harbor Seal, Herring, Chinook",
                   "2.Harbor Seal,  Herring, Hake", "3.Harbor Seal,  Herring, Smolts", "4.Harbor Seal,  Chinook, Hake",
                   "6. Harbor Seal, Smolts, Hake",
                   "Int","Int Sm")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}

model.selectionPREY.1 <- model.selection.PREY(dataPrey.1, n, dataPrey.1$TP)
x<-data.frame(model.selectionPREY.1)
subset(x, delAICc<=2)
modelPREY1<-lmer(TP~Location.2+HakeBiomass+(1|AA), data=dataPrey.1)
summary(modelPREY1)

model.selectionPREY.2 <- model.selection.PREY(dataPrey.2, n, dataPrey.2$TP)
x<-data.frame(model.selectionPREY.2)
subset(x, delAICc<=2)
modelPREY2<-lmer(TP~Location.2+allSmolt+(1|AA), data=dataPrey.2)
summary(modelPREY2)

model.selectionPREY.3 <- model.selection.PREY(dataPrey.3, n, dataPrey.3$TP)
x<-data.frame(model.selectionPREY.3)
subset(x, delAICc<=2)
modelPREY3<-lmer(TP~Location.2+allSmolt+HakeBiomass+Herring.Biomass+(1|AA), data=dataPrey.3)
summary(modelPREY3)

model.selectionPREY.4 <- model.selection.PREY(dataPrey.4, n, dataPrey.4$TP)
x<-data.frame(model.selectionPREY.4)
subset(x, delAICc<=5)
modelPREY3<-lmer(TP~Location.2+allSmolt+HakeBiomass+Herring.Biomass+(1|AA), data=dataPrey.3)
summary(modelPREY3)

model.selectionPREY.5 <- model.selection.PREY(dataPrey.5, n, dataPrey.5$TP)
x<-data.frame(model.selectionPREY.5)
subset(x, delAICc<=5)
modelPREY3<-lmer(TP~Location.2+allSmolt+HakeBiomass+Herring.Biomass+(1|AA), data=dataPrey.3)
summary(modelPREY3)

########################    Environmental  Plots#######################
library(dotwhisker)
library(broom)
library(dplyr)
library(colorspace)
library(ggpubr)

color<- rep(c('black','#CCA65A','#7EBA68','#6FB1E7','#D494E1'), each=4)
color2<- rep(c('black','#CCA65A','#7EBA68','#6FB1E7','#D494E1'), times=4)
color3<- rep(c('black','#CCA65A','#7EBA68','#6FB1E7','#D494E1'), each=3)
color4<- rep(c('black','#CCA65A','#7EBA68','#6FB1E7','#D494E1'), times=3)
color5<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=2)
color6<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=2)
color7<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=5)
color8<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=5)



ENV <- tidy(modelENV.1)
x <- data.frame("term" = c('Intercept',   'Location','Summer Upwelling'), "estimate" =   c(ENV[1:3,2]), "std.error" = c(ENV[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
ENV.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.1)$AA[1,1]),
                "std.error" = c(ranef(modelENV.1)$AA[1,1]), "group" = 'random', model="Ala")
ENV.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.1)$AA[3,1]) ,
                "std.error" = c(ranef(modelENV.1)$AA[3,1]), "group" = 'random', model="Pro")
ENV.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.1)$AA[2,1]), 
                "std.error" = c(ranef(modelENV.1)$AA[2,1]), "group" = 'random', model="Glu")
ENV.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.1)$AA[4,1]), 
                "std.error" = c(ranef(modelENV.1)$AA[4,1]), "group" = 'random', model="Val")
ENV.VAL <- as_tibble(x)

ENV.mods <- rbind(ENV.fixed, ENV.GLU, ENV.ALA,   ENV.VAL, ENV.PRO)

ENV.plot1 <-small_multiple(ENV.mods) +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("A. Year-0") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10)) 
ENV.plot1

ENV <- tidy(modelENV.2)
x <- data.frame("term" = c('Intercept',   'Location','Summer SST'), "estimate" =   c(ENV[1:3,2]), "std.error" = c(ENV[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
ENV.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.2)$AA[1,1]),
                "std.error" = c(ranef(modelENV.2)$AA[1,1]), "group" = 'random', model="Ala")
ENV.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.2)$AA[3,1]) ,
                "std.error" = c(ranef(modelENV.2)$AA[3,1]), "group" = 'random', model="Pro")
ENV.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.2)$AA[2,1]), 
                "std.error" = c(ranef(modelENV.2)$AA[2,1]), "group" = 'random', model="Glu")
ENV.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.2)$AA[4,1]), 
                "std.error" = c(ranef(modelENV.2)$AA[4,1]), "group" = 'random', model="Val")
ENV.VAL <- as_tibble(x)

ENV.mods <- rbind(ENV.fixed, ENV.GLU, ENV.ALA,   ENV.VAL, ENV.PRO)

ENV.plot2 <-small_multiple(ENV.mods) +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("B. Year-1") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10)) 
ENV.plot2


ENV <- tidy(modelENV.3)
x <- data.frame("term" = c('Intercept',   'Location','Summer Upwelling', "MEI"), "estimate" =   c(ENV[1:4,2]), 
                "std.error" = c(ENV[1:4,3]), "group" = c(rep('fixed', 4)), model="fixed")
ENV.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.3)$AA[1,1]),
                "std.error" = c(ranef(modelENV.3)$AA[1,1]), "group" = 'random', model="Ala")
ENV.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.3)$AA[3,1]) ,
                "std.error" = c(ranef(modelENV.3)$AA[3,1]), "group" = 'random', model="Pro")
ENV.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.3)$AA[2,1]), 
                "std.error" = c(ranef(modelENV.3)$AA[2,1]), "group" = 'random', model="Glu")
ENV.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV.3)$AA[4,1]), 
                "std.error" = c(ranef(modelENV.3)$AA[4,1]), "group" = 'random', model="Val")
ENV.VAL <- as_tibble(x)

ENV.mods <- rbind(ENV.fixed, ENV.GLU, ENV.ALA,   ENV.VAL, ENV.PRO)

ENV.plot3 <-small_multiple(ENV.mods) +
  theme_bw()+
  geom_point (colour=color)+
  geom_errorbar (colour=color2, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("C. Year-2") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10)) 
ENV.plot3


ENV.plot<-ggarrange(ENV.plot1,ENV.plot2,ENV.plot3,
                     ncol = 3, nrow = 1, align= 'v')
ENV.plot<-annotate_figure(ENV.plot,
                           top = text_grob("Ocean Condition Models", color = "black", face = "bold", size = 14)
)
ENV.plot


pdf(file="Results/Figures/HCoefPlotENV.pdf", width=12, height=6.5)
ENV.plot

dev.off()

########################    Foodweb Plots#######################


PREY <- tidy(modelPREY1)
x <- data.frame("term" = c('Intercept',  'Location', "Hake"), "estimate" = fixef(modelPREY1), "std.error" = c(PREY[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
PREY.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY1)$AA[1,1]), "std.error" = c(ranef(modelPREY1)$AA[1,1]),   "group" = 'random', model="Ala")
PREY.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY1)$AA[2,1]), "std.error" = c(ranef(modelPREY1)$AA[2,1]), "group" = 'random', model="Glu")
PREY.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY1)$AA[3,1]), "std.error" = c(ranef(modelPREY1)$AA[3,1]), "group" = 'random', model="Pro")
PREY.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY1)$AA[4,1]), "std.error" = c(ranef(modelPREY1)$AA[4,1]), "group" = 'random', model="Val")
PREY.VAL <- as_tibble(x)
PREY.mods <- rbind(PREY.fixed, PREY.GLU, PREY.ALA,  PREY.VAL, PREY.PRO)


PREY.plot1 <-small_multiple(PREY.mods) +
  theme_bw( )+
 # theme(plot.margin = margin(1, 0.35, 0.5, 0.35, "cm"))+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("A. Year-0") +
  # scale_y_discrete("", breaks=waiver(), labels=waiver())+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))
PREY.plot1

PREY <- tidy(modelPREY2)
x <- data.frame("term" = c('Intercept',  'Location', "Chinook Smolts"), "estimate" = fixef(modelPREY2), "std.error" = c(PREY[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
PREY.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY2)$AA[1,1]), "std.error" = c(ranef(modelPREY2)$AA[1,1]),   "group" = 'random', model="Ala")
PREY.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY2)$AA[2,1]), "std.error" = c(ranef(modelPREY2)$AA[2,1]), "group" = 'random', model="Glu")
PREY.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY2)$AA[3,1]), "std.error" = c(ranef(modelPREY2)$AA[3,1]), "group" = 'random', model="Pro")
PREY.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY2)$AA[4,1]), "std.error" = c(ranef(modelPREY2)$AA[4,1]), "group" = 'random', model="Val")
PREY.VAL <- as_tibble(x)
PREY.mods <- rbind(PREY.fixed, PREY.GLU, PREY.ALA,  PREY.VAL, PREY.PRO)


PREY.plot2 <-small_multiple(PREY.mods) +
  theme_bw()+
  geom_point (colour=color3)+
  #theme(plot.margin = margin(1, 0.35, 0.5, 0.35, "cm"))+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("B. Year-1") +
  # scale_y_discrete("", breaks=waiver(), labels=waiver())+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))
PREY.plot2

PREY <- tidy(modelPREY3)
x <- data.frame("term" = c('Intercept',  'Location', "Chinook Smolts", "Hake", "Herring"), "estimate" = fixef(modelPREY3), "std.error" = c(PREY[1:5,3]), "group" = c(rep('fixed', 5)), model="fixed")
PREY.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY3)$AA[1,1]), "std.error" = c(ranef(modelPREY3)$AA[1,1]),   "group" = 'random', model="Ala")
PREY.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY3)$AA[2,1]), "std.error" = c(ranef(modelPREY3)$AA[2,1]), "group" = 'random', model="Glu")
PREY.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY3)$AA[3,1]), "std.error" = c(ranef(modelPREY3)$AA[3,1]), "group" = 'random', model="Pro")
PREY.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY3)$AA[4,1]), "std.error" = c(ranef(modelPREY3)$AA[4,1]), "group" = 'random', model="Val")
PREY.VAL <- as_tibble(x)
PREY.mods <- rbind(PREY.fixed, PREY.GLU, PREY.ALA,  PREY.VAL, PREY.PRO)


PREY.plot3 <-small_multiple(PREY.mods) +
  theme_bw()+
  geom_point (colour=color7)+
  geom_errorbar (colour=color8, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
 # theme(plot.margin = margin(1, 0.35, 0.5, 0.35, "cm"))+
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("C. Year-2") +
  # scale_y_discrete("", breaks=waiver(), labels=waiver())+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))
PREY.plot3

PREY.plot<-ggarrange(PREY.plot1,PREY.plot2,PREY.plot3,
          ncol = 3, nrow = 1, align= 'v')
PREY.plot<-annotate_figure(PREY.plot,
                top = text_grob("Food Web Models", color = "black", face = "bold", size = 14)
)
PREY.plot


pdf(file="Results/Figures/HCoefPlot.pdf", width=12, height=6.5)
PREY.plot
dev.off()





pdf(file="Results/Figures/HCoefPlot.pdf", width=5, height=8)
ggarrange(ENV.plot, PREY.plot,
          ncol = 1, nrow = 2, align= 'v')

dev.off()




###############Residual Plots ##############
library(HLMdiag)
HLMresid(modelFULL, level=1)
col<- c('#CCA65A','#7EBA68','#6FB1E7','#D494E1')



dataresid<-cbind(dataCLIM,resid=resid(modelENV))


dataresid.coastal<- subset(dataresid, Location.2=="Coastal")
dataresid.coastal<- subset(dataresid.coastal, AA=="VAL"|AA=="GLU"|AA=="ALA"|AA=="PRO")

#dataresid.coastal <- dataresid.coastal[complete.cases(dataresid.coastal), ]
dataresid.coastal$AA<- factor(dataresid.coastal$AA, levels = c("GLU","ALA","VAL",
                                                               "PRO"),
                              labels = c("Glutamic Acid", "Alanine", "Valine", "Proline")
)


resid.coastalENV <- ggplot(dataresid.coastal, aes(x = Year, y = resid)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  labs(y="Residuals")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=4),aes(color = AA, alpha=0.5))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(dataresid.coastal$AA))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)
resid.coastal






dataresid.ss<- subset(dataresid, Location.2=="Inland")
dataresid.ss<- subset(dataresid.ss, AA=="VAL"|AA=="GLU"|AA=="ALA"|AA=="PRO")

dataresid.ss <- dataresid.ss[complete.cases(dataresid.ss), ]
dataresid.ss$AA<- factor(dataresid.ss$AA, levels = c("GLU","ALA",
                                                     "VAL","PRO"),
                         labels = c("Glutamic Acid", "Alanine", "Valine", "Proline")
)


resid.ssENV <- ggplot(dataresid.ss, aes(x = Year, y = resid)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  labs(y="Residuals")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  geom_smooth(method="gam", aes(color = AA, alpha=0.5))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(dataresid.ss$AA))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)
resid.ss


dataresid<-cbind(dataPrey,resid=resid(modelPREY))




dataresid.coastal<- subset(dataresid, Location.2=="Coastal")
dataresid.coastal<- subset(dataresid.coastal, AA=="VAL"|AA=="GLU"|AA=="ALA"|AA=="PRO")

#dataresid.coastal <- dataresid.coastal[complete.cases(dataresid.coastal), ]
dataresid.coastal$AA<- factor(dataresid.coastal$AA, levels = c("GLU","ALA","VAL",
                                                               "PRO"),
                              labels = c("Glutamic Acid", "Alanine", "Valine", "Proline")
)


resid.coastalPREY <- ggplot(dataresid.coastal, aes(x = Year, y = resid)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  labs(y="Residuals")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k=4),aes(color = AA, alpha=0.5))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(dataresid.coastal$AA))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)
resid.coastal






dataresid.ss<- subset(dataresid, Location.2=="Inland")
dataresid.ss<- subset(dataresid.ss, AA=="VAL"|AA=="GLU"|AA=="ALA"|AA=="PRO")

dataresid.ss <- dataresid.ss[complete.cases(dataresid.ss), ]
dataresid.ss$AA<- factor(dataresid.ss$AA, levels = c("GLU","ALA",
                                                     "VAL","PRO"),
                         labels = c("Glutamic Acid", "Alanine", "Valine", "Proline")
)


resid.ssPREY <- ggplot(dataresid.ss, aes(x = Year, y = resid)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  labs(y="Residuals")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  geom_smooth(method="gam", aes(color = AA, alpha=0.5))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(dataresid.ss$AA))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)
resid.ss

dev.off()

