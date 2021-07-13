rm(list = ls())

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

########################    Hier Clim2 Models       ############################

dataCLIM <- subset(data.lag1, AA=="GLU")
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

dataCLIM <- subset(data.lag2, AA=="GLU")
dataCLIM.2 <-dataCLIM.2 %>% select(MEI, 
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

dataCLIM<- subset(data.lag3, AA=="GLU")
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

dataCLIM<- subset(data.lag4, AA=="GLU")
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
  
  aic.output <- rbind(
                      AICc(lm(y~Location.2, data=dataframe)), #2
                      AICc(lm(y~PDO+Location.2, data=dataframe)),#3
                      AICc(lm(y~NPGO+Location.2, data=dataframe)),#4
                      AICc(lm(y~MEI+Location.2, data=dataframe)),#5
                      AICc(lm(y~UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      AICc(lm(y~PDO+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#7
                      #AICc(lm(y~PDO+MEI+Location.2, data=dataframe)),
                      AICc(lm(y~UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#8
                      #AICc(lm(y~MEI+NPGO+Location.2, data=dataframe)),#5
                      AICc(lm(y~MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#9
                      #AICc(lm(y~MEI+NPGO+UpInAn.45.Spring+Location.2, data=dataframe)),#7
                      #AICc(lm(y~PDO+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#8
                      #AICc(lm(y~PDO+MEI+NPGO+Location.2, data=dataframe)),#9
                      #AICc(lm(y~PDO+NPGO+UpInAn.45.Spring+Location.2, data=dataframe)),
                      
                      AICc(lm(y~Location.2+WA.SST.Su, data=dataframe)), #10
                      AICc(lm(y~WA.SST.Su+PDO+Location.2, data=dataframe)),#11
                      AICc(lm(y~WA.SST.Su+NPGO+Location.2, data=dataframe)),#12
                      AICc(lm(y~WA.SST.Su+MEI+Location.2, data=dataframe)),#13
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Spring+Location.2, data=dataframe)),#14
                      # AIC(lm(y~WA.SST.Su+PDO+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#15
                      #AIC(lm(y~WA.SST.Su+PDO+MEI+Location.2, data=dataframe)),
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#16
                      #AIC(lm(y~WA.SST.Su+MEI+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#17
                      AICc(lm(y~Location.2+UpInAn.45.Summer, data=dataframe)), #18
                      AICc(lm(y~UpInAn.45.Summer+PDO+Location.2, data=dataframe)),#19
                      AICc(lm(y~UpInAn.45.Summer+NPGO+Location.2, data=dataframe)),#20
                      AICc(lm(y~UpInAn.45.Summer+MEI+Location.2, data=dataframe)),#21
                      #AIC(lm(y~UpInAn.45.Summer+UpInAn.45.Spring+Location.2, data=dataframe)),#22
                      # AIC(lm(y~UpInAn.45.Summer+PDO+NPGO+Location.2, data=dataframe)),
                      #AICc(lm(y~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),
                      #AIC(lm(y~UpInAn.45.Summer+PDO+MEI+Location.2, data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#23
                      #AIC(lm(y~UpInAn.45.Summer+MEI+NPGO+Location.2, data=dataframe)),
                      #AICc(lm(y~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#24
                      
                      AICc(lm(y~Location.2+Col.Dis.high, data=dataframe)), #25
                      AICc(lm(y~Col.Dis.high+PDO+Location.2, data=dataframe)),#26
                      AICc(lm(y~Col.Dis.high+NPGO+Location.2, data=dataframe)),#27
                      AICc(lm(y~Col.Dis.high+MEI+Location.2, data=dataframe)),#28
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Spring+Location.2, data=dataframe)),#29
                      # AICc(lm(y~Col.Dis.high+PDO+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#30
                      #AICc(lm(y~Col.Dis.high+PDO+MEI+Location.2, data=dataframe)),
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#31
                      # AICc(lm(y~Col.Dis.high+MEI+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#32
                      
                      AICc(lm(y~Col.Dis.high+WA.SST.Su+Location.2, data=dataframe)),#33
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Summer+Location.2, data=dataframe)),#34
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Summer+Location.2, data=dataframe)),#35
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2, data=dataframe)),
                      AICc(lm(y~MEI*Location.2, data=dataframe))#36
                      #36
                      
                      
                      
  )
  
  model.names <- c("2. Location", "3. PDO", "4. NPGO", "5. MEI", "6. Upwelling (Sp)","NPGO,PDO", "7. PDO, Upwelling (Sp)",  
                   "8. NPGO, Upwelling (Sp)","9. MEI, Upwelling (Sp)",  
                   
                   "10. WA.SST.Su, Location", "11. WA.SST.Su, PDO", "12. WA.SST.Su, NPGO", "13. WA.SST.Su, MEI",
                   "14. WA.SST.Su, Upwelling (Sp)", "15. WA.SST.Su, PDO, Upwelling (Sp)", 
                   "16. WA.SST.Su, NPGO, Upwelling (Sp)", "17. WA.SST.Su, MEI, Upwelling (Sp)",
                   
                   "18. UpInAn.45.Summer, Location", "19. UpInAn.45.Summer, PDO", "20. UpInAn.45.Summer, NPGO", 
                   "21. UpInAn.45.Summer, MEI", #"22. UpInAn.45.Summer, Upwelling (Sp)",
                   #"23.UpInAn.45.Summer,  PDO, Upwelling (Sp)",
                   "24. UpInAn.45.Summer,  NPGO, Upwelling (Sp)", #"25. UpInAn.45.Summer, MEI, Upwelling (Sp)",
                   
                   "Col.Dis.high, Location", "Col.Dis.high, PDO", "Col.Dis.high, NPGO", "Col.Dis.high, MEI", "Upwelling (Sp)",
                   "2.Col.Dis.high,  PDO, Upwelling (Sp)", "4.Col.Dis.high,  NPGO, Upwelling (Sp)",
                   "6. Col.Dis.high, MEI, Upwelling (Sp)",
                   
                   "Col.Dis.high, WA.SST.Su", " Col.Dis.high UpInAn.45.Summer", "WA.SST.Su, UpInAn.45.Summer", "WA.SST.Su, UpInAn.45.Summer, Col.Dis.high",
                   "PDO INT")
  
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
modelENV.1 <- lm(TP~Location.2+UpInAn.45.Summer, data=dataCLIM.1)
summary(modelENV.1)

model.selectionCLIM.2 <- ModelSelection.CLIM(dataCLIM.2, n, dataCLIM.2$TP)
clim.2 <-data.frame(model.selectionCLIM.2)
subset(clim.2, delAICc<=1.97)
modelENV.2 <- lm(TP~Location.2+WA.SST.Su, data=dataCLIM.2)
summary(modelENV.2)

model.selectionCLIM.3 <- ModelSelection.CLIM(dataCLIM.3, n, dataCLIM.3$TP)
clim.3 <-data.frame(model.selectionCLIM.3)
subset(clim.3, delAICc<=2)
modelENV.3 <- lm(TP~Location.2+UpInAn.45.Summer+MEI, data=dataCLIM.3)
summary(modelENV.3)

model.selectionCLIM.4 <- ModelSelection.CLIM(dataCLIM.4, n, dataCLIM.4$TP)
clim.4 <-data.frame(model.selectionCLIM.4)
subset(clim.4, delAICc<=2)
modelENV.4 <- lm(TP~Location.2+WA.SST.Su, data=dataCLIM.4)
summary(modelENV.4)

########################    Hier PREY Models       ############################

dataPrey <- subset(data.lag1, AA=="GLU")
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

dataPrey <- subset(data.lag2, AA=="GLU")
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

dataPrey <- subset(data.lag3, AA=="GLU")
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

dataPrey <- subset(data.lag4, AA=="GLU")
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

dataPrey <- subset(data.lag4, AA=="GLU")
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
  
  aic.output <- rbind( 
                      AICc(lm(y~Location.2, data=dataframe)),
                      AICc(lm(y~Herring.Biomass+Location.2, data=dataframe)),
                      AICc(lm(y~Chinook+Location.2, data=dataframe)),
                      AICc(lm(y~allSmolt+Location.2, data=dataframe)),
                      AICc(lm(y~HakeBiomass+Location.2, data=dataframe)),
                      AICc(lm(y~Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AICc(lm(y~Herring.Biomass+HakeBiomass, data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AICc(lm(y~HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AICc(lm(y~allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AICc(lm(y~allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      AICc(lm(y~allSmolt+Chinook+HakeBiomass+Location.2, data=dataframe)),#7
                      AICc(lm(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2, data=dataframe)),#8
                      AICc(lm(y~Herring.Biomass+allSmolt+Chinook+Location.2, data=dataframe)),#9
                      AICc(lm(y~Herring.Biomass+Chinook+HakeBiomass+Location.2, data=dataframe)),
                      
                      
                      AICc(lm(y~HarborSeal+Location.2, data=dataframe)), 
                      AICc(lm(y~HarborSeal+Herring.Biomass+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+Chinook+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+allSmolt+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+HakeBiomass+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AICc(lm(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AICc(lm(y~HarborSeal+Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AICc(lm(y~HarborSeal+HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      #AICc(lm(y~HarborSeal+allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AICc(lm(y~HarborSeal+allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      
                      
                      AICc(lm(y~Herring.Biomass*HakeBiomass+Location.2, data=dataframe))

                      
                      
  )
  
  model.names <- c("Location", "Herring", "Chinook", "Smolts", "Hake","1. Herring, Chinook", 
                   "2. Herring, Hake", "3. Herring, Smolts", "4. Chinook, Hake","CH SM",
                   "6. Smolts, Hake", "CH SM", "8. Herring, Smolts, Hake","CH SM",
                   "11. Herring Hake Chinook",
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook", "Harbor Seal, Smolts", 
                   "Hake","1. Harbor Seal, Herring, Chinook",
                   "2.Harbor Seal,  Herring, Hake", "3.Harbor Seal,  Herring, Smolts", "4.Harbor Seal,  Chinook, Hake",
                   "6. Harbor Seal, Smolts, Hake",
                   "Int")
  
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
modelPREY1<-lm(TP~Location.2+HakeBiomass, data=dataPrey.1)
summary(modelPREY1)

model.selectionPREY.2 <- model.selection.PREY(dataPrey.2, n, dataPrey.2$TP)
x<-data.frame(model.selectionPREY.2)
subset(x, delAICc<=2)
modelPREY2<-lm(TP~Location.2+allSmolt, data=dataPrey.2)
summary(modelPREY2)

model.selectionPREY.3 <- model.selection.PREY(dataPrey.3, n, dataPrey.3$TP)
x<-data.frame(model.selectionPREY.3)
subset(x, delAICc<=2)
modelPREY3<-lm(TP~Location.2+allSmolt+HakeBiomass+Herring.Biomass, data=dataPrey.3)
summary(modelPREY3)

model.selectionPREY.4 <- model.selection.PREY(dataPrey.4, n, dataPrey.4$TP)
x<-data.frame(model.selectionPREY.4)
subset(x, delAICc<=5)
modelPREY3<-lm(TP~Location.2+allSmolt+HakeBiomass+Herring.Biomass, data=dataPrey.3)
summary(modelPREY3)

model.selectionPREY.5 <- model.selection.PREY(dataPrey.5, n, dataPrey.5$TP)
x<-data.frame(model.selectionPREY.5)
subset(x, delAICc<=5)
modelPREY3<-lm(TP~Location.2+allSmolt+HakeBiomass+Herring.Biomass, data=dataPrey.3)
summary(modelPREY3)
