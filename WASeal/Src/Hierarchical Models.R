
library(mgcViz)
library(AICcmodavg)
library(vapniks)
require(nlme)
require(mgcv)
library(qpcR)
library(dplyr)
library(dotwhisker)
library(car)
library(ggpubr)

library(devtools)
#install_version("lme4", "1.1-25")
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


dataALL <-  read.csv('Data/Compiled/HierarchicalData1AdCorr.csv')
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag1ad<-data

dataALL <-  read.csv('Data/Compiled/HierarchicalData2AdCorr.csv')
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag2ad<-data

dataALL <-  read.csv('Data/Compiled/HierarchicalData3AdCorr.csv')
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag3ad<-data

dataALL <-  read.csv('Data/Compiled/HierarchicalData0.csv')
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data.lag3ad<-data
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

dataCLIM<- subset(data.lag1ad, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataCLIM.1ad <-dataCLIM %>% select(MEI, 
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
dataCLIM.1ad <- dataCLIM.1ad[complete.cases(dataCLIM.1ad), ]

dataCLIM<- subset(data.lag2ad, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataCLIM.2ad <-dataCLIM %>% select(MEI, 
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
dataCLIM.2ad <- dataCLIM.2ad[complete.cases(dataCLIM.2ad), ]

dataCLIM<- subset(data.lag3ad, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataCLIM.3ad <-dataCLIM %>% select(MEI, 
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
dataCLIM.3ad <- dataCLIM.3ad[complete.cases(dataCLIM.3ad), ]


dataCLIM<- subset(data.lag0, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataCLIM.0 <-dataCLIM %>% select(MEI, 
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
dataCLIM.0 <- dataCLIM.0[complete.cases(dataCLIM.0), ]

ModelSelection.CLIM<- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe)), #1
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe)), #2
                      AICc(lmer(y~PDO+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~MEI+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      AICc(lmer(y~PDO+NPGO+Location.2+(1|AA), data=dataframe)),#7
                      AICc(lmer(y~PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#8
                      AICc(lmer(y~UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#9
                      AICc(lmer(y~MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#10
                      AICc(lmer(y~Location.2+(1|AA)+WA.SST.Su, data=dataframe)), #11
                      AICc(lmer(y~WA.SST.Su+PDO+Location.2+(1|AA), data=dataframe)),#12
                      AICc(lmer(y~WA.SST.Su+NPGO+Location.2+(1|AA), data=dataframe)),#13
                      AICc(lmer(y~WA.SST.Su+MEI+Location.2+(1|AA), data=dataframe)),#14
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#15
                      AICc(lmer(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#16
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#17
                      AICc(lmer(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#18
                      AICc(lmer(y~Location.2+(1|AA)+UpInAn.45.Summer, data=dataframe)), #19
                      AICc(lmer(y~UpInAn.45.Summer+PDO+Location.2+(1|AA), data=dataframe)),#20
                      AICc(lmer(y~UpInAn.45.Summer+NPGO+Location.2+(1|AA), data=dataframe)),#21
                      AICc(lmer(y~UpInAn.45.Summer+MEI+Location.2+(1|AA), data=dataframe)),#22
                      AICc(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#23
                      AICc(lmer(y~Location.2+(1|AA)+Col.Dis.high, data=dataframe)), #24
                      AICc(lmer(y~Col.Dis.high+PDO+Location.2+(1|AA), data=dataframe)),#25
                      AICc(lmer(y~Col.Dis.high+NPGO+Location.2+(1|AA), data=dataframe)),#26
                      AICc(lmer(y~Col.Dis.high+MEI+Location.2+(1|AA), data=dataframe)),#27
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#28
                      AICc(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#29
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#30
                      AICc(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#31
                      AICc(lmer(y~Col.Dis.high+WA.SST.Su+Location.2+(1|AA), data=dataframe)),#32
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#33
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#34
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2+(1|AA), data=dataframe))#35
                      
                      
                      
                      
  )
  
  model.names <- c("1. Null",
                   "2. Location", 
                   "3. PDO", 
                   "4. NPGO", 
                   "5. MEI", 
                   "6. Upwelling (Spring)",
                   "7. NPGO,PDO", 
                   "8. PDO, Upwelling (Spring)",  
                   "9. NPGO, Upwelling (Spring)",
                   "10. MEI, Upwelling (Spring)",  
                   "11. SST (Summer)", 
                   "12. SST (Summer), PDO", 
                   "13. SST (Summer), NPGO",
                   "14. SST (Summer), MEI",
                   "15. SST (Summer), Upwelling (Spring)",
                   "16. SST (Summer), PDO, Upwelling (Spring)", 
                   "17. SST (Summer), NPGO, Upwelling (Spring)",
                   "18. SST (Summer), MEI, Upwelling (Spring)",
                   
                   "19. Upwelling (Summer)", 
                   "20. Upwelling (Summer), PDO", 
                   "21. Upwelling (Summer), NPGO", 
                   "22. Upwelling (Summer), MEI",
                   "23. Upwelling (Summer),  NPGO, Upwelling (Spring)", 
                   "24. Columbia Discharge High, Location", 
                   "25. Columbia Discharge High, PDO", 
                   "26. Columbia Discharge High, NPGO",
                   "27. Columbia Discharge High, MEI", 
                   "28. Upwelling (Spring)",
                   "29.Columbia Discharge High,  PDO, Upwelling (Spring)",
                   "30. Columbia Discharge High,  NPGO, Upwelling (Spring)",
                   "31. Columbia Discharge High, MEI, Upwelling (Spring)",
                   
                   "32. Columbia Discharge High, SST (Summer)", 
                   "33. Columbia Discharge High, Upwelling (Summer)", 
                   "34. SST (Summer), Upwelling (Summer)", 
                   "35. SST (Summer), Upwelling (Summer), Columbia Discharge High")
  
  # row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- round(aic.weight1/sum(aic.weight1), digits=2)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(model.names, round(aic.output, digits=2), round(delaic, digits=2), aic.weight)
  
  colnames(aic.output)<- c("Covariates","AICc", "delAICc", "AICc Weight")
  return(aic.output)
}

model.selectionCLIM.0 <- ModelSelection.CLIM(dataCLIM.0, n, dataCLIM.0$TP)
clim.0 <-data.frame(model.selectionCLIM.0)
clim.0.ordered <- clim.1[order(clim.0$AICc),]
subset(clim.0, delAICc<=5)
clim.1.ordered[1:6,]

model.selectionCLIM.1 <- ModelSelection.CLIM(dataCLIM.1, n, dataCLIM.1$TP)
clim.1 <-data.frame(model.selectionCLIM.1)
clim.1.ordered <- clim.1[order(clim.1$AICc),]
subset(clim.1, delAICc<=5)
clim.1.ordered[1:5,]

sjPlot::tab_df(clim.1.ordered[1:5,],
               title = "Physiological Delay Top 5 Models", 
               file = "Results/Tables/Clim1Top5.doc")

modelENV.1 <- lmer(TP~Location.2+UpInAn.45.Summer+(1|AA), data=dataCLIM.1)
summary(modelENV.1)
modelENV.1 <- lmer(TP~Location.2+UpInAn.45.Summer+(1|AA), data=dataCLIM.1)
summary(lmer(TP~(1|AA), data=dataCLIM.1))

model.selectionCLIM.2 <- ModelSelection.CLIM(dataCLIM.2, n, dataCLIM.2$TP)
clim.2 <-data.frame(model.selectionCLIM.2)
clim.2.ordered <- clim.2[order(clim.2$AICc),]
subset(clim.2, delAICc<=1.97)
clim.2.ordered[1:5,]

sjPlot::tab_df(clim.2.ordered[1:5,],
               title = "1-Year Ecological Delay Top 5 Models", 
               file = "Results/Tables/Clim2Top5.doc")


modelENV.2 <- lmer(TP~Location.2+WA.SST.Su+(1|AA), data=dataCLIM.2)
summary(modelENV.2)

model.selectionCLIM.3 <- ModelSelection.CLIM(dataCLIM.3, n, dataCLIM.3$TP)
clim.3 <-data.frame(model.selectionCLIM.3)
subset(clim.3, delAICc<=2)
clim.3.ordered <- clim.3[order(clim.3$AICc),]
clim.3.ordered[1:6,]

sjPlot::tab_df(clim.3.ordered[1:5,],
               title = "1-Year Ecological Delay Top 5 Models", 
               file = "Results/Tables/Clim3Top5.doc")



modelENV.3 <- lmer(TP~Location.2+Col.Dis.high+(1|AA), data=dataCLIM.3)
summary(modelENV.3)

model.selectionCLIM.4 <- ModelSelection.CLIM(dataCLIM.4, n, dataCLIM.4$TP)
clim.4 <-data.frame(model.selectionCLIM.4)
subset(clim.4, delAICc<=2)
modelENV.4 <- lmer(TP~Location.2+WA.SST.Su+(1|AA), data=dataCLIM.4)
summary(modelENV.4)

model.selectionCLIM.1ad <- ModelSelection.CLIM(dataCLIM.1ad, n, dataCLIM.1ad$TP)
clim.1ad <-data.frame(model.selectionCLIM.1ad)
subset(clim.1ad, delAICc<=2)
clim.1ad.ordered <- clim.1ad[order(clim.1ad$AICc),]
clim.1ad.ordered[1:6,]

model.selectionCLIM.2ad <- ModelSelection.CLIM(dataCLIM.2ad, n, dataCLIM.2ad$TP)
clim.2ad <-data.frame(model.selectionCLIM.2ad)
subset(clim.2ad, delAICc<=2)
clim.2ad.ordered <- clim.2ad[order(clim.2ad$AICc),]
clim.2ad.ordered[1:6,]

model.selectionCLIM.3ad <- ModelSelection.CLIM(dataCLIM.3ad, n, dataCLIM.3ad$TP)
clim.3ad <-data.frame(model.selectionCLIM.3ad)
subset(clim.3ad, delAICc<=2)
clim.3ad.ordered <- clim.3ad[order(clim.3ad$AICc),]
clim.3ad.ordered[1:6,]
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

dataPREY<- subset(data.lag1ad, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataPREY.1ad <-dataPREY %>% select(allSmolt, 
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
dataPREY.1ad <- dataPREY.1ad[complete.cases(dataPREY.1ad), ]

dataPREY<- subset(data.lag2ad, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataPREY.2ad <-dataPREY %>% select(allSmolt, 
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
dataPREY.2ad <- dataPREY.2ad[complete.cases(dataPREY.2ad), ]

dataPREY<- subset(data.lag3ad, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataPREY.3ad <-dataPREY %>% select(allSmolt, 
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
dataPREY.3ad <- dataPREY.3ad[complete.cases(dataPREY.3ad), ]

dataPREY<- subset(data.0, AA=="GLU"|AA=="ALA"|AA=="PRO"|AA=="VAL"|AA=="PRO")
dataPREY.0 <-dataPREY %>% select(allSmolt, 
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
dataPREY.0 <- dataPREY.0[complete.cases(dataPREY.0), ]

model.selection.PREY <- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe)), #1, null
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe)), #2 Location Only
                      AICc(lmer(y~Location.2+WildProduction+(1|AA), data=dataframe)), #3. Wild Production Smolts
                      #AICc(lmer(y~Location.2+HatcherySmolts+(1|AA), data=dataframe)),#4. Hatchery smolts
                      AICc(lmer(y~Herring.Biomass+Location.2+(1|AA), data=dataframe)), #5. Herring Biomass
                      AICc(lmer(y~Chinook+Location.2+(1|AA), data=dataframe)),#6. Chinook Escapements
                      AICc(lmer(y~allSmolt+Location.2+(1|AA), data=dataframe)), #7. All smolts
                      AICc(lmer(y~HakeBiomass+Location.2+(1|AA), data=dataframe)), #8. Hake Biomass
                      AICc(lmer(y~Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#9. Herring Biomass, Chinook Escapements
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+(1|AA), data=dataframe)),#10. Herring Biomass, Hake Biomass
                      AICc(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#11. Herring Biomass, all smolt
                      AICc(lmer(y~HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#12. Hake Bimass, Chinook Escapement
                      AICc(lmer(y~allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#13. all smolt, Chinook Escapements
                      AICc(lmer(y~allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#14. all smolt, Hake Biomass
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),#15. all smolt, Chinook escapent, Hake Biomas
                      AICc(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#16. Herring Biomass, all smolts, Hake Biomass
                      AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#17.Herring biomass, all smolts, Chinook escapements
                      AICc(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)), #18. Herring Biomass, Chinook escapements, Hake biomass
                      AICc(lmer(y~HarborSeal+Location.2+(1|AA), data=dataframe)), #19. Harbor seal 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA), data=dataframe)),#20. Harbor seal, herring biomass
                      AICc(lmer(y~HarborSeal+Chinook+Location.2+(1|AA), data=dataframe)), #21. Harbor seal, Chinook escapement
                      AICc(lmer(y~HarborSeal+allSmolt+Location.2+(1|AA), data=dataframe)),#22. Harbor seal, all smolt
                      AICc(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|AA), data=dataframe)),#23. Harbor seal, hake biomass
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#24. Harbor seal, herring biomass, Chinook escapement
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#25. Harbor seal, herring biomass, hakr biomass
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#26. Harbor seal, herring biomass, all smolts
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#27. Harbor seal, hake biomass, Chinool escapement
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#28. Harbor seal, all smolts, hake biomass
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)), #29. Hake:Herring biomass
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+allSmolt+Location.2+(1|AA), data=dataframe))#30. Hake:Herring biomass, all smolt
                      
                      
                      
  )
  
  model.names <- c("1. Null",
                   "2. Location Only",
                   "3. Herring Biomass", 
                   "4. Chinook Escapements", 
                   "5. All Chinook Smolts", 
                   "6. Hake Biomass",
                   "7. Herring Biomass, Chinook escapements", 
                   "8. Herring Biomass, Hake Biomass", 
                   "9. Herring Biomass, All Chinook Smolts", 
                   "10. Chinook escapement, Hake biomass",
                   "11. Chinook escapments, All Chinook Smolts",
                   "12. All Chinook Smolts, Hake biomass", 
                   "13. Chinook escapement, All Chinook Smolts, Hake biomass", 
                   "14. Herring Biomass, All Chinook Smolts, Hake biomass",
                   "15. Chinook escapement, All Chinook Smolts, herring biomass",
                   "16. Herring biomass, Hake biomass, Chinook escapement",
                   "17. Harbor Seal population",
                   "18. Harbor Seal population, Herring biomass", 
                   "19. Harbor Seal population, Chinook escapement", 
                   "20.Harbor Seal population, All Chinook Smolts", 
                   "21. Harbor seal population, Hake biomass",
                   "22. Harbor Seal population, Herring biomass, Chinook escapement",
                   "23.Harbor Seal population,  Herring biomass, Hake biomass", 
                   "24.Harbor Seal population,  Herring biomass, All Chinook Smolts", 
                   "25.Harbor Seal population,  Chinook escapement, Hake biomass",
                   "26. Harbor Seal population, All Chinook Smolts, Hake biomass",
                   "27. Hake and Herring interaction",
                   "28. Hake and Herring interaction, All Chinook Smolts")
  
  #row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(model.names, round(aic.output, digits=2), round(delaic, digits=2), round(aic.weight, digits=2))
  
  colnames(aic.output)<- c("Covariates","AICc", "delAICc", "AICc Weight")
  return(aic.output)
}


model.selectionPREY.0 <- model.selection.PREY(dataPREY.0, n, dataPREY.0$TP)
prey0<-data.frame(model.selectionPREY.0)
subset(prey0, delAICc<=2)
prey0.ordered <- prey0[order(prey0$AICc),]
prey0.ordered[1:6,]

model.selectionPREY.1 <- model.selection.PREY(dataPrey.1, n, dataPrey.1$TP)
prey1<-data.frame(model.selectionPREY.1)
subset(prey1, delAICc<=2)
prey1.ordered <- prey1[order(prey1$AICc),]
prey1.ordered[1:6,]
modelPREY1<-lmer(TP~Location.2+Chinook+Herring.Biomass+(1|AA), data=dataPrey.1)
summary(modelPREY1)
sjPlot::tab_df(prey1.ordered[1:5,],
               title = "Physiological Delay Top 5 Models (Prey Availability)", 
               file = "Results/Tables/Prey1Top5.doc")


model.selectionPREY.2 <- model.selection.PREY(dataPrey.2, n, dataPrey.2$TP)
prey2<-data.frame(model.selectionPREY.2)
prey2.ordered <- prey2[order(prey2$AICc),]
prey2.ordered[1:6,]
modelPREY2<-lmer(TP~Location.2+HakeBiomass+(1|AA), data=dataPrey.2)
summary(modelPREY2)
sjPlot::tab_df(prey2.ordered[1:5,],
               title = "1-Year Physiological Top 5 Models (Prey Availability)", 
               file = "Results/Tables/Prey2Top5.doc")



model.selectionPREY.3 <- model.selection.PREY(dataPrey.3, n, dataPrey.3$TP)
prey3<-data.frame(model.selectionPREY.3)
prey3.ordered <- prey3[order(prey3$AICc),]
prey3.ordered[1:7,]
modelPREY3<-lmer(TP~Location.2+allSmolt+(1|AA), data=dataPrey.3)
summary(modelPREY3)
sjPlot::tab_df(prey3.ordered[1:7,],
               title = "2-Year Physiological Top 5 Models (Prey Availability)", 
               file = "Results/Tables/Prey3Top5.doc")




model.selectionPREY.1ad <- model.selection.PREY(dataPREY.1ad, n, dataPREY.1ad$TP)
prey1ad<-data.frame(model.selectionPREY.1ad)
subset(prey1ad, delAICc<=5)
subset(prey1ad, delAICc<=2)
prey1ad.ordered <- prey1ad[order(prey1ad$AICc),]
prey1ad.ordered[1:6,]



model.selectionPREY.2ad <- model.selection.PREY(dataPREY.2ad, n, dataPREY.2ad$TP)
prey2ad<-data.frame(model.selectionPREY.2ad)
prey2ad.ordered <- prey2ad[order(prey2ad$AICc),]
prey2ad.ordered[1:6,]
summary(lmer(TP~Location.2+HarborSeal+Herring.Biomass+(1|AA), dataPREY.2ad))
summary(lmer(TP~Location.2+HarborSeal+(1|AA), dataPREY.2ad))
summary(lmer(TP~Location.2+HarborSeal+HakeBiomass+(1|AA), dataPREY.2ad))


sjPlot::tab_df(prey2ad.ordered[1:5,],
               title = "2-Year Future Top 5 Models (Prey Availability)", 
               file = "Results/Tables/FutureTop5.doc")



model.selectionPREY.3ad <- model.selection.PREY(dataPREY.3ad, n, dataPREY.3ad$TP)
prey3ad<-data.frame(model.selectionPREY.3ad)
prey3ad.ordered <- prey3ad[order(prey3ad$AICc),]
prey3ad.ordered[1:6,]


######################## EDITED plots for VERION 2 #######################
color.env <-c('#7EBA68','#6FB1E7','#D494E1',"#009ADE")
plot_summs(modelENV.1, modelENV.2, modelENV.3, 
           model.names = c("Physiological Delay", "1-year Ecological Delay", "2-year Ecological Delay"), colors=color.env, ylab="Covariates")+
  theme_bw()+
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Ocean Condition Models") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        
        #axis.text.x = element_blank(),
        axis.text.y = element_text( size = 10)) 



plot_summs(modelPREY1, modelPREY2, modelPREY3, 
           model.names = c("Physiological Delay", "1-year Ecological Delay", "2-year Ecological Delay"), colors=color.env, ylab="Covariates")+
  theme_bw()+geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web Models") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        
        #axis.text.x = element_blank(),
        axis.text.y = element_text( size = 10)) 

########################    1. Environmental  Plots#######################
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
x <- data.frame("term" = c('Intercept',   'Location','Summer Upwelling'), "estimate" =   c(ENV[1:3,4]), "std.error" = c(ENV[1:3,5]), "group" = c(rep('fixed', 3)), model="fixed")
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

ENV.mods <- rbind(ENV.fixed, ENV.GLU, ENV.ALA,ENV.VAL, ENV.PRO)

ENV.plot1 <-dwplot(ENV.mods) +
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
x <- data.frame("term" = c('Intercept',   'Location','Summer SST'), "estimate" =   c(ENV[1:3,4]), "std.error" = c(ENV[1:3,5]), "group" = c(rep('fixed', 3)), model="fixed")
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
x <- data.frame("term" = c('Intercept',   'Location','Summer Upwelling', "MEI"), "estimate" =   c(ENV[1:4,4]), 
                "std.error" = c(ENV[1:4,5]), "group" = c(rep('fixed', 4)), model="fixed")
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


########################    1. Foodweb Plots#######################


PREY <- tidy(modelPREY1)
x <- data.frame("term" = c('Intercept',  'Location', "Hake"), "estimate" = fixef(modelPREY1), "std.error" = c(PREY[1:3,5]), "group" = c(rep('fixed', 3)), model="fixed")
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


PREY.plot<-ggarrange(PREY.plot3,PREY.plot2,PREY.plot1,
                     ncol = 3, nrow = 1, align= 'v')
PREY.plot<-annotate_figure(PREY.plot,
                           top = text_grob("Food Web Models", color = "black", face = "bold", size = 14)
)
PREY.plot




########################    2. Environmental  Plots#######################
library(dotwhisker)
library(broom)
library(dplyr)
library(colorspace)
library(ggpubr)
#'black','#CCA65A','#7EBA68','#6FB1E7','#D494E1',"#009ADE",

colore1bar<- c('black','grey','grey',"grey",'grey', #Int
               '#CCA65A','grey','grey',"grey",'grey', #Location
               '#D494E1','grey','grey',"grey",'grey') #Upwelling

colore1pt<- c( 'black','#CCA65A','#D494E1','grey',"grey",'grey', 'grey',
               'grey', 'grey','grey','grey',"grey",'grey', 'grey',
               'grey') #Location


ENV <- tidy(modelENV.1)
x <- data.frame("term" = c('Intercept',   'Location','Upwelling'), "estimate" =   c(ENV[1:3,4]), "std.error" = c(ENV[1:3,5]), "group" = c(rep('fixed', 3)), model="fixed")
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
  geom_point (colour=colore1pt)+
  geom_errorbar (colour=colore1bar, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Physiological Delay") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text( size = 10)) 
ENV.plot1

ENV <- tidy(modelENV.2)
x <- data.frame("term" = c('Intercept',   'Location','Summer SST'), "estimate" =   c(ENV[1:3,4]), "std.error" = c(ENV[1:3,5]), "group" = c(rep('fixed', 3)), model="fixed")
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



colore2bar<- c('black','grey','grey',"grey",'grey', #Int
               '#CCA65A','grey','grey',"grey",'grey', #Location
               '#6FB1E7','grey','grey',"grey",'grey') #SST

colore2pt<- c( 'black','#CCA65A','#6FB1E7','grey',"grey",'grey', 'grey',
               'grey', 'grey','grey','grey',"grey",'grey', 'grey',
               'grey') #Location


ENV.plot2 <-small_multiple(ENV.mods) +
  theme_bw()+
  geom_point (colour=colore2pt)+
  geom_errorbar (colour=colore2bar, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("1-year Ecological Delay") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text( size = 10)) 
ENV.plot2


ENV <- tidy(modelENV.3)
x <- data.frame("term" = c('Intercept',   'Location','Upwelling', "MEI"), "estimate" =   c(ENV[1:4,4]), 
                "std.error" = c(ENV[1:4,5]), "group" = c(rep('fixed', 4)), model="fixed")
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


colore3bar<- c('black','grey','grey',"grey",'grey', #Int
               '#CCA65A','grey','grey',"grey",'grey', #Location
               '#7EBA68','grey','grey',"grey",'grey',#MEI
               '#D494E1','grey','grey',"grey",'grey') #upwelling
colore3pt<- c( 'black','#CCA65A','#7EBA68','#D494E1',
               'grey',"grey",'grey', 'grey',
               'grey', 'grey','grey','grey',
               "grey",'grey', 'grey','grey',
               'grey','grey','grey','grey') #Location


ENV.plot3 <-small_multiple(ENV.mods) +
  theme_bw()+
  geom_point (colour=colore3pt)+
  geom_errorbar (colour=colore3bar, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("2-year Ecological Delay") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text( size = 10)) 
ENV.plot3


ENV.plot<-ggarrange(ENV.plot1,ENV.plot2,ENV.plot3,
                    ncol = 3, nrow = 1, align= 'v')
ENV.plot<-annotate_figure(ENV.plot,
                          top = text_grob("Ocean Condition Models", color = "black", face = "bold", size = 14)
)
ENV.plot


pdf(file="Results/Figures/HCoefPlot_nonweightedBeta.rev.pdf", width=12, height=6.5)
ENV.plot
dev.off()

########################   2. Foodweb Plots#######################
colorp1bar<- c('#7EBA68','grey','grey',"grey",'grey', #Hake
               'black','grey','grey',"grey",'grey', #Int
               '#CCA65A','grey','grey',"grey",'grey') #Location

colorp1pt<- c('#7EBA68', 'black','#CCA65A','grey',"grey",'grey', 'grey',
              'grey', 'grey','grey','grey',"grey",'grey', 'grey',
              'grey') #Location

PREY <- tidy(modelPREY1)
x <- data.frame("term" = c('Intercept',  'Location', "Hake"), "estimate" = fixef(modelPREY1), "std.error" = c(PREY[1:3,5]), "group" = c(rep('fixed', 3)), model="fixed")
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
  geom_point (colour=colorp1pt)+
  geom_errorbar (colour=colorp1bar, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle(" ") +
  # scale_y_discrete("", breaks=waiver(), labels=waiver())+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))
PREY.plot1



PREY <- tidy(modelPREY2)
x <- data.frame("term" = c('Intercept',  'Location', "Smolts"), "estimate" = fixef(modelPREY2), "std.error" = c(PREY[1:3,5]), "group" = c(rep('fixed', 3)), model="fixed")
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

colorp2bar<- c('black','grey','grey',"grey",'grey', #chinook
               '#CCA65A','grey','grey',"grey",'grey', #Int
               '#6FB1E7','grey','grey',"grey",'grey') #Location

colorp2pt<- c( 'black','#CCA65A','#6FB1E7','grey',"grey",'grey', 'grey',
               'grey', 'grey','grey','grey',"grey",'grey', 'grey',
               'grey') #Location

PREY.plot2 <-small_multiple(PREY.mods) +
  theme_bw()+
  geom_point (colour=colorp2pt)+
  #theme(plot.margin = margin(1, 0.35, 0.5, 0.35, "cm"))+
  geom_errorbar (colour=colorp2bar, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web Models") +
  # scale_y_discrete("", breaks=waiver(), labels=waiver())+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=16), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))
PREY.plot2

PREY <- tidy(modelPREY3)
x <- data.frame("term" = c('Intercept',  'Location', "Smolts", "Hake", "Herring"), "estimate" = fixef(modelPREY3), "std.error" = c(PREY[1:5,5]), "group" = c(rep('fixed', 5)), model="fixed")
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

labs<- c("Chinook Smolts","Hake","Herring",'Intercept','Location')

colorp3bar<- c('#7EBA68','grey','grey',"grey",'grey', #chinook
               '#D494E1','grey','grey',"grey",'grey', #Hake
               'black' ,'grey','grey',"grey",'grey', #herring
               '#CCA65A','grey','grey',"grey",'grey', #Int
               '#6FB1E7','grey','grey',"grey",'grey') #Location


colorp3pt<- c('#7EBA68', '#D494E1','black','#CCA65A','#6FB1E7', 
              'grey',"grey",'grey', 'grey', 'grey',
              'grey', 'grey','grey','grey',"grey",
              'grey', 'grey','grey','grey',"grey",
              'grey','grey','grey','grey',"grey") #Location


PREY.plot3 <-small_multiple(PREY.mods) +
  theme_bw()+
  geom_point (colour=colorp3pt)+
  geom_errorbar (colour=colorp3bar, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle(" ") +
  # facet_grid(term ~ ., labeller(term=label_wrap_gen(10)))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))
PREY.plot3

PREY.plot<-ggarrange(PREY.plot1,PREY.plot2,PREY.plot3,
                     ncol = 3, nrow = 1, align= 'v')
PREY.plot<-annotate_figure(PREY.plot
                           # top = text_grob("Food Web Models", color = "black", face = "bold", size = 14)
)
PREY.plot


pdf(file="Results/Figures/Prey_HCoefPlot_nonweightedBeta.rev.pdf", width=12, height=6.5)
PREY.plot
dev.off()


####FULL PLOT####


FULL.plot<-ggarrange(ENV.plot3,ENV.plot2,ENV.plot1,
                     PREY.plot3,PREY.plot2,PREY.plot1,
                     ncol = 3, nrow = 2, align= 'v')
FULL.plot<-annotate_figure(FULL.plot,
                           top = text_grob("Ocean Condition Models", color = "black", face = "bold", size = 16)
)
FULL.plot


pdf(file="Results/Figures/HCoefPlot.FULL.pdf", width=9.5, height=8.5)
FULL.plot
dev.off()



###############Residual Plots ##############
library(HLMdiag)
HLMresid(modelFULL, level=1)
col<-c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7')
palette(c('#7EBA68','#00C1B2','#CCA65A','#6FB1E7','#D494E1'))

pdf(file="Results/Figures/ResidualsDiagnosticsENV.pdf", width=16, height=12)
par(mfrow=c(3,4))

y<- dataCLIM.1
plot(predict(modelENV.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="A. Physiological Delay Model", col=dataCLIM.1$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelENV.1), residuals(modelENV.1), ylab = "Residuals", xlab = "Predicted",
     main= "Physiological Delay Model",  pch=19, cex=0.5, col=dataCLIM.1$AA)
abline(0, 0)

y<- dataCLIM.1
plot(predict(modelENV.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="B. Physiological Delay Model", col=dataCLIM.1$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelENV.1), residuals(modelENV.1), ylab = "Residuals", xlab = "Predicted",
     main= "Physiological Delay Model",  pch=19, cex=0.5, col=dataCLIM.1$AA)
abline(0, 0)

y<- dataCLIM.1
plot(predict(modelENV.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="C. Physiological Delay Model", col=dataCLIM.1$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelENV.1), residuals(modelENV.1), ylab = "Residuals", xlab = "Predicted",
     main= "Physiological Delay Model",  pch=19, cex=0.5, col=dataCLIM.1$AA)
abline(0, 0)

dev.off()
#######TABLES######

tab_df(prey1.ordered,
       
       title = "Prey Physiological Delay", #always give
       #your tables
       #titles
       file = "Results/Tables/Prey1.doc")

tab_df(prey2.ordered,
       
       title = "1-year Ecological Delay", #always give
       #your tables
       #titles
       file = "Results/Tables/Prey2.doc")
tab_df(prey3.ordered,
       
       title = "1-year Ecological Delay", #always give
       #your tables
       #titles
       file = "Results/Tables/Prey3.doc")

tab_df(clim.1.ordered,
       
       title = "Environment Physiological Delay", #always give
       #your tables
       #titles
       file = "Results/Tables/Env1.doc")

tab_df(clim.2.ordered,
       
       title = "1-year Environmental Delay", #always give
       #your tables
       #titles
       file = "Results/Tables/Env2.doc")
tab_df(clim.3.ordered,
       
       title = "2-year Environmental Delay", #always give
       #your tables
       #titles
       file = "Results/Tables/Env3.doc")
