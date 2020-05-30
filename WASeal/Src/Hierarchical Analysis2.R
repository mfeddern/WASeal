library(lme4)
library(AICcmodavg)
data <-  read.csv("Data/Compiled/HierarchicalData2.csv")

########################    Hier Clim Models       ############################

dataCLIM <- subset(data, Year>=1950&Year<=2008 & AA=="Glu"|AA=="PRO"|AA=="ALA")
#dataCLIM <- subset(data, Year>=1960&Year<=2008)# & AA=="Glu"|AA=="PRO"|AA=="ALA")

dataCLIM <-dataCLIM %>% select(MEI, 
                            PDO,
                            NPGO,
                            WA.SST.Su, 
                            UpInAn.45.Spring,
                            Col.Dis.high,
                            UpInAn.45.Summer,
                            TP.norm,
                            TP,
                            Year,
                            AA,
                            Sample.ID,
                            Location.2)



dataCLIM <- dataCLIM[complete.cases(dataCLIM), ]
ModelSelection.CLIM<- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~PDO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~PDO+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~PDO+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~PDO+MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~UpInAn.45.Spring+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~MEI+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~MEI+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AIC(lmer(y~MEI+NPGO+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#7
                      AIC(lmer(y~PDO+MEI+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#8
                      AIC(lmer(y~PDO+MEI+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#9
                      AIC(lmer(y~PDO+NPGO+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      
                      AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID)+WA.SST.Su, data=dataframe)), 
                      AIC(lmer(y~WA.SST.Su+PDO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+PDO+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~WA.SST.Su+PDO+MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~WA.SST.Su+MEI+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID)+UpInAn.45.Summer, data=dataframe)), 
                      AIC(lmer(y~UpInAn.45.Summer+PDO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+PDO+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~UpInAn.45.Summer+PDO+MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~UpInAn.45.Summer+MEI+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID)+Col.Dis.high, data=dataframe)), 
                      AIC(lmer(y~Col.Dis.high+PDO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+PDO+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~Col.Dis.high+PDO+MEI+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~Col.Dis.high+MEI+NPGO+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~Col.Dis.high+WA.SST.Su+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~Col.Dis.high+UpInAn.45.Summer+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Summer+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2+(1|AA)+(1|Sample.ID), data=dataframe))#6
                      
                      
                      
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
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}
n<- 17
model.selectionCLIM <- ModelSelection.CLIM(dataCLIM, n, dataCLIM$TP)
model.selectionCLIM
length(model.selectionCLIM[,1])
x<-data.frame(model.selectionCLIM)
subset(x, delAICc<=1.97)
modelENV<- lmer(TP~Location.2+Col.Dis.high+(1|AA)+(1|Sample.ID), data=dataCLIM)
summary(modelENV)




ModelSelection.CLIM2<- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lmer(y~Location.2+(1|AA), data=dataframe)), 
                      AIC(lmer(y~PDO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~NPGO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~MEI+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~PDO+NPGO+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~PDO+MEI+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~MEI+NPGO+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(y~MEI+NPGO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#7
                      AIC(lmer(y~PDO+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#8
                      AIC(lmer(y~PDO+MEI+NPGO+Location.2+(1|AA), data=dataframe)),#9
                      AIC(lmer(y~PDO+NPGO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      
                      AIC(lmer(y~Location.2+(1|AA)+WA.SST.Su, data=dataframe)), 
                      AIC(lmer(y~WA.SST.Su+PDO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+NPGO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+MEI+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~WA.SST.Su+PDO+NPGO+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~WA.SST.Su+PDO+MEI+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~WA.SST.Su+MEI+NPGO+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+(1|AA)+UpInAn.45.Summer, data=dataframe)), 
                      AIC(lmer(y~UpInAn.45.Summer+PDO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+NPGO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+MEI+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~UpInAn.45.Summer+PDO+NPGO+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~UpInAn.45.Summer+PDO+MEI+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~UpInAn.45.Summer+MEI+NPGO+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+(1|AA)+Col.Dis.high, data=dataframe)), 
                      AIC(lmer(y~Col.Dis.high+PDO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+NPGO+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+MEI+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Col.Dis.high+PDO+NPGO+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~Col.Dis.high+PDO+MEI+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~Col.Dis.high+MEI+NPGO+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      
                      AIC(lmer(y~Col.Dis.high+WA.SST.Su+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~Col.Dis.high+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2+(1|AA), data=dataframe))#6
                      
                      
                      
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
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}
n<- 17
model.selectionCLIM2 <- ModelSelection.CLIM2(dataCLIM, n, dataCLIM$TP.norm)
model.selectionCLIM2
x<-data.frame(model.selectionCLIM)
modelENV<- lmer(TP~Location.2+Col.Dis.high+(1|AA)+(1|Sample.ID), data=dataCLIM)
summary(modelENV)


########################    Hier PREY Models       ############################

data2 <- subset(data, Year>=1973&Year<=2008&AA=="Glu"|AA=="PRO"|AA=="ALA")
#data2 <- subset(data, Year>=1973&Year<=2008)#&AA=="Glu"|AA=="PRO"|AA=="ALA")
dataPrey <-data2 %>% select(allSmolt, 
                            HakeBiomass,
                            Herring.Biomass,
                            AA,
                            Chinook, 
                            HarborSeal,
                            Coho,
                            Chum,
                            TP.norm,
                            TP,
                            Year,
                            Sample.ID,
                            Location.2)
dataPrey <- dataPrey[complete.cases(dataPrey), ]
ModelSelection.WAPREY <- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lmer(y~(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~Herring.Biomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      # AIC(lmer(y~allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Herring.Biomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      # AIC(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      #  AIC(lmer(y~allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      #  AIC(lmer(y~allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      #  AIC(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#7
                      # AIC(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#8
                      #AIC(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#9
                      AIC(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      
                      AIC(lmer(y~Location.2+Chum+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~Chum+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chum+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      # AIC(lmer(y~Chum+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chum+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chum+Herring.Biomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Chum+Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      # AIC(lmer(y~Chum+Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~Chum+HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      #AIC(lmer(y~Chum+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      #AIC(lmer(y~Chum+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      # AIC(lmer(y~Location.2+Coho+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~Coho+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Coho+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      #AIC(lmer(y~Coho+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Coho+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Coho+Herring.Biomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Coho+Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      #AIC(lmer(y~Coho+Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~Coho+HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      # AIC(lmer(y~Coho+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      # AIC(lmer(y~Coho+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+HarborSeal+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HarborSeal+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      # AIC(lmer(y~HarborSeal+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      # AIC(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      # AIC(lmer(y~HarborSeal+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      # AIC(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~HarborSeal+Chum+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~HarborSeal+Coho+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~Chum+Coho+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AIC(lmer(y~Chum+Coho+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe))#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  
  
  model.names <- c("Location", "Herring", "Chinook", "Hake","1. Herring, Chinook", "2. Herring, Hake",  "4. Chinook, Hake",
                   
                   "11. Herring Hake Chinook",
                   
                   "Chum, Location", "Chum, Herring", "Chum, Chinook",  "Chum, Hake","1. Chum,  Herring, Chinook", "2. Chum, Herring, Hake", 
                   "4. Chum, Chinook, Hake",
                   "Coho, Herring", "Coho, Chinook",  "Coho, Hake","1. Coho, Herring, Chinook", "2.Coho,  Herring, Hake",
                   "4.Coho,  Chinook, Hake", 
                   
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook",  "Hake","1. Harbor Seal, Herring, Chinook",
                   "2.Harbor Seal,  Herring, Hake", "4.Harbor Seal,  Chinook, Hake", 
                   
                   
                   "Harbor Seal, Chum", " Harbor Seal Coho", "Chum, Coho", "Chum, Coho, Harbor Seal")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}
n<- 16
ModelSelection.WAPREY <- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AICc(lmer(y~Herring.Biomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Herring.Biomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass+allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AICc(lmer(y~HakeBiomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AICc(lmer(y~allSmolt+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AICc(lmer(y~allSmolt+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#7
                      AICc(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#8
                      AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#9
                      AICc(lmer(y~Herring.Biomass+Chinook+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      
                      AICc(lmer(y~(1|Location.2)+Chum+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AICc(lmer(y~Chum+Herring.Biomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Chum+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Chum+allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Chum+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Chum+Herring.Biomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AICc(lmer(y~Chum+Herring.Biomass+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AICc(lmer(y~Chum+Herring.Biomass+allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AICc(lmer(y~Chum+HakeBiomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AICc(lmer(y~Chum+allSmolt+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AICc(lmer(y~Chum+allSmolt+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AICc(lmer(y~(1|Location.2)+Coho+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AICc(lmer(y~Coho+Herring.Biomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Coho+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Coho+allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Coho+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~Coho+Herring.Biomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AICc(lmer(y~Coho+Herring.Biomass+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AICc(lmer(y~Coho+Herring.Biomass+allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AICc(lmer(y~Coho+HakeBiomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AICc(lmer(y~Coho+allSmolt+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AICc(lmer(y~Coho+allSmolt+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AICc(lmer(y~(1|Location.2)+HarborSeal+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~HarborSeal+allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~HarborSeal+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+allSmolt+Chinook+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AICc(lmer(y~HarborSeal+Chum+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+Coho+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AICc(lmer(y~Chum+Coho+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AICc(lmer(y~Chum+Coho+HarborSeal+(1|Location.2)+(1|AA)+(1|Sample.ID), data=dataframe))#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
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
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}

ModelSelection.WAPREY <- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~Location.2+(1|Sample.ID/AA), data=dataframe)), 
                      AICc(lmer(y~Herring.Biomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Herring.Biomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#1
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),#3
                      AICc(lmer(y~HakeBiomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#4
                      AICc(lmer(y~allSmolt+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#5
                      AICc(lmer(y~allSmolt+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#6
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#7
                      AICc(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#8
                      AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#9
                      AICc(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      
                      AICc(lmer(y~Location.2+Chum+(1|Sample.ID/AA), data=dataframe)), 
                      AICc(lmer(y~Chum+Herring.Biomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Chum+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Chum+allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Chum+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Chum+Herring.Biomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#1
                      AICc(lmer(y~Chum+Herring.Biomass+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#2
                      AICc(lmer(y~Chum+Herring.Biomass+allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),#3
                      AICc(lmer(y~Chum+HakeBiomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#4
                      AICc(lmer(y~Chum+allSmolt+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#5
                      AICc(lmer(y~Chum+allSmolt+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#6
                      
                      AICc(lmer(y~Location.2+Coho+(1|Sample.ID/AA), data=dataframe)), 
                      AICc(lmer(y~Coho+Herring.Biomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Coho+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Coho+allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Coho+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~Coho+Herring.Biomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#1
                      AICc(lmer(y~Coho+Herring.Biomass+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#2
                      AICc(lmer(y~Coho+Herring.Biomass+allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),#3
                      AICc(lmer(y~Coho+HakeBiomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#4
                      AICc(lmer(y~Coho+allSmolt+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#5
                      AICc(lmer(y~Coho+allSmolt+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#6
                      
                      AICc(lmer(y~Location.2+HarborSeal+(1|Sample.ID/AA), data=dataframe)), 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#1
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#2
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|Sample.ID/AA), data=dataframe)),#3
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+allSmolt+Chinook+Location.2+(1|Sample.ID/AA), data=dataframe)),#5
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|Sample.ID/AA), data=dataframe)),#6
                      
                      AICc(lmer(y~HarborSeal+Chum+Location.2+(1|Sample.ID/AA), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+Coho+Location.2+(1|Sample.ID/AA), data=dataframe)),#5
                      AICc(lmer(y~Chum+Coho+Location.2+(1|Sample.ID/AA), data=dataframe)),#6
                      AICc(lmer(y~Chum+Coho+HarborSeal+Location.2+(1|Sample.ID/AA), data=dataframe))#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
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
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}


model.selectionPREY <- ModelSelection.WAPREY(dataPrey, n, dataPrey$TP)
length(model.selectionPREY[,1])
model.selectionPREY
x<-data.frame(model.selectionPREY)
subset(x, delAICc<=1.97)
summary(lmer(TP~Location.2+Herring.Biomass+(1|Sample.ID/AA), data=dataPrey))
lmer(TP~Location.2+Herring.Biomass+(1|AA), data=dataPrey)



ModelSelection.WAPREY2 <- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lmer(y~Location.2+(1|AA), data=dataframe)), 
                      AIC(lmer(y~Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Chinook+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~allSmolt+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),#7
                      AIC(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#8
                      AIC(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#9
                      AIC(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      
                      AIC(lmer(y~Location.2+Chum+(1|AA), data=dataframe)), 
                      AIC(lmer(y~Chum+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Chum+Chinook+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Chum+allSmolt+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Chum+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Chum+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~Chum+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~Chum+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~Chum+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~Chum+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~Chum+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+Coho+(1|AA), data=dataframe)), 
                      AIC(lmer(y~Coho+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Coho+Chinook+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Coho+allSmolt+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Coho+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Coho+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~Coho+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~Coho+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~Coho+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~Coho+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~Coho+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+HarborSeal+(1|AA), data=dataframe)), 
                      AIC(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~HarborSeal+Chinook+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~HarborSeal+allSmolt+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~HarborSeal+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      
                      AIC(lmer(y~HarborSeal+Chum+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~HarborSeal+Coho+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~Chum+Coho+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(y~Chum+Coho+HarborSeal+Location.2+(1|AA), data=dataframe))#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Herring", "Chinook",  "Hake","1. Herring, Chinook", "2. Herring, Hake", "4. Chinook, Hake",
                    
                   "11. Herring Hake Chinook",
                   
                   "Chum, Location", "Chum, Herring", "Chum, Chinook",  "Chum, Hake","1. Chum,  Herring, Chinook", "2. Chum, Herring, Hake", 
                   "4. Chum, Chinook, Hake",
                   
                   "Coho, Location", "Coho, Herring",  "Coho, Hake","1. Coho, Herring, Chinook", "2.Coho,  Herring, Hake",
                   "4.Coho,  Chinook, Hake", 
                   
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook",  "Hake","1. Harbor Seal, Herring, Chinook",
                   "2.Harbor Seal,  Herring, Hake", "4.Harbor Seal,  Chinook, Hake",
                  
                   
                   "Harbor Seal, Chum", " Harbor Seal Coho", "Chum, Coho", "Chum, Coho, Harbor Seal")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}
n<- 16

model.selectionPREY <- ModelSelection.WAPREY2(dataPrey, n, dataPrey$TP)
length(model.selectionPREY[,1])
model.selectionPREY
x<-data.frame(model.selectionPREY)
subset(x, delAICc<=1.97)

summary(lmer(TP~HarborSeal+Location.2+(1|AA), data=dataPrey))

########################     Hierarchical Nutrient Models       ############################
dataNut <- subset(data, Year>=1928&Year<=2014 & AA=="Glu"|AA=="PRO"|AA=="ALA")
dataNut <-dataNut %>% select(PHE.mean,
                            PHE.norm,
                            AA,
                            d13C.s, 
                            d13C.norm,
                            TP.norm,
                            TP,
                            Year,
                            Sample.ID,
                            Location.2)
dataNut <- dataNut[complete.cases(dataNut), ]



ModelSelection.WANut <- function(dataframe,n, y) {
  
  aic.output <- rbind(
    AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
    AIC(lmer(y~PHE.norm+d13C.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
    AIC(lmer(y~PHE.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
    AIC(lmer(y~d13C.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe))
  )
  
  names <- seq(1,n,1)
  model.names <- c( "Location","Phe, 13C", "PHE", "13C")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}

n<-4
model.selectionNUTRIENT <- ModelSelection.WANut(dataNut, n, dataNut$TP.norm)
model.selectionNUTRIENT
summary(lmer(TP.norm~PHE.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataNut))








########################     Plots#######################
library(dotwhisker)
library(broom)
library(dplyr)
library(colorspace)
library(ggpubr)


ENV <- tidy(modelENV)
x <- data.frame("term" = c('Intercept',  'Climate', 'Location'), "estimate" = fixef(modelENV), "std.error" = c(ENV[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
ENV.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[1,1],coef(modelENV)$AA[1,3]) , "std.error" = c(ranef(modelENV)$AA[1,1], ranef(modelENV)$AA[1,2]), "group" = 'random', model="Glu")
ENV.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[2,1],coef(modelENV)$AA[2,3]), "std.error" = c(ranef(modelENV)$AA[2,1], ranef(modelENV)$AA[1,2]), "group" = 'random', model="Asp")
ENV.ASP <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[3,1],coef(modelENV)$AA[3,3]), "std.error" = c(ranef(modelENV)$AA[3,1],  ranef(modelENV)$AA[1,2]), "group" = 'random', model="Pro")
ENV.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[4,1],coef(modelENV)$AA[4,3]), "std.error" = c(ranef(modelENV)$AA[4,1], ranef(modelENV)$AA[1,2]), "group" = 'random', model="Val")
ENV.VAL <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[5,1],coef(modelENV)$AA[5,3]), "std.error" = c(ranef(modelENV)$AA[5,1],  ranef(modelENV)$AA[1,2]), "group" = 'random', model="Ala")
ENV.ALA <- as_tibble(x)
ENV.mods <- rbind(ENV.fixed, ENV.GLU, ENV.ASP, ENV.PRO, ENV.VAL, ENV.ALA)



NUTRIENT <- tidy(modelNUTRIENT)
x <- data.frame("term" = c('Intercept','Phenylalanine', 'd13C','Location'), "estimate" = fixef(modelNUTRIENT),
                "std.error" = c(NUTRIENT[1:4,3]), "group" = c(rep('fixed', 4)), model="fixed")
NUTRIENT.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept', "Phenylalanine", 'Location'), "estimate" = c(coef(modelNUTRIENT)$AA[1,1],coef(modelNUTRIENT)$AA[1,2],coef(modelNUTRIENT)$AA[1,4]), 
                "std.error" = c(ranef(modelNUTRIENT)$AA[1,1], ranef(modelNUTRIENT)$AA[1,3], ranef(modelNUTRIENT)$AA[1,2]), "group" = 'random', model="Glu")
NUTRIENT.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept', "Phenylalanine", 'Location'), "estimate" = c(coef(modelNUTRIENT)$AA[2,1],coef(modelNUTRIENT)$AA[2,2],coef(modelNUTRIENT)$AA[2,3]),
                "std.error" = c(ranef(modelNUTRIENT)$AA[2,1], ranef(modelNUTRIENT)$AA[2,3], ranef(modelNUTRIENT)$AA[2,2]), "group" = 'random', model="Asp")
NUTRIENT.ASP <- as_tibble(x)
x <- data.frame("term" = c('Intercept', "Phenylalanine", 'Location'), "estimate" = c(coef(modelNUTRIENT)$AA[3,1],coef(modelNUTRIENT)$AA[3,2],coef(modelNUTRIENT)$AA[3,3]), 
                "std.error" = c(ranef(modelNUTRIENT)$AA[3,1],  ranef(modelNUTRIENT)$AA[3,3],  ranef(modelNUTRIENT)$AA[3,2]), "group" = 'random', model="Pro")
NUTRIENT.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept', "Phenylalanine", 'Location'), "estimate" = c(coef(modelNUTRIENT)$AA[4,1],coef(modelNUTRIENT)$AA[4,2],coef(modelNUTRIENT)$AA[4,3]), 
                "std.error" = c(ranef(modelNUTRIENT)$AA[4,1], ranef(modelNUTRIENT)$AA[4,3], ranef(modelNUTRIENT)$AA[4,2]), "group" = 'random', model="Val")
NUTRIENT.VAL <- as_tibble(x)
x <- data.frame("term" = c('Intercept', "Phenylalanine", 'Location'), "estimate" = c(coef(modelNUTRIENT)$AA[5,1],coef(modelNUTRIENT)$AA[5,2],coef(modelNUTRIENT)$AA[5,3]),
                "std.error" = c(ranef(modelNUTRIENT)$AA[5,1],  ranef(modelNUTRIENT)$AA[5,3],ranef(modelNUTRIENT)$AA[5,2]), "group" = 'random', model="Ala")
NUTRIENT.ALA <- as_tibble(x)
NUTRIENT.mods <- rbind(NUTRIENT.fixed, NUTRIENT.GLU, NUTRIENT.ASP, NUTRIENT.PRO, NUTRIENT.VAL, NUTRIENT.ALA)


modelPREY <- lmer(TP~Herring.Biomass+HarborSeal+Location.2+(1|AA), data=data3)
PREY <- tidy(modelPREY)
x <- data.frame("term" = c('Intercept',  'Herring', 'Harbor Seal', 'Location'), "estimate" = fixef(modelPREY), 
                "std.error" = c(PREY[1:4,3]), "group" = c(rep('fixed', 4)), model="fixed")
PREY.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[1,1]) ,
                "std.error" = c(ranef(modelPREY)$AA[1,1]), "group" = 'random', model="Glu")
PREY.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[2,1]),
                "std.error" = c(ranef(modelPREY)$AA[2,1]), "group" = 'random', model="Asp")
PREY.ASP <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[3,1]), 
                "std.error" = c(ranef(modelPREY)$AA[3,1]), "group" = 'random', model="Pro")
PREY.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[4,1]), 
                "std.error" = c(ranef(modelPREY)$AA[4,1]), "group" = 'random', model="Val")
PREY.VAL <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[5,1]), 
                "std.error" = c(ranef(modelPREY)$AA[5,1]), "group" = 'random', model="Ala")
PREY.ALA <- as_tibble(x)
PREY.mods <- rbind(PREY.fixed, PREY.GLU, PREY.ASP, PREY.PRO, PREY.VAL, PREY.ALA)





FULL <- tidy(modelFULL)
x <- data.frame("term" = c('Intercept', "Pheylalanine", "d13C", 'Harbor Seal', 'Location'), "estimate" = fixef(modelFULL), 
                "std.error" = c(FULL[1:5,3]), "group" = c(rep('fixed', 5)), model="fixed")
FULL.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Pheylalanine',  "Location"), "estimate" = c(coef(modelFULL)$AA[1,1],coef(modelFULL)$AA[1,2],coef(modelFULL)$AA[1,5]) ,
                "std.error" = c(ranef(modelFULL)$AA[1,1], ranef(modelFULL)$AA[1,3], ranef(modelFULL)$AA[1,2]), "group" = 'random', model="Glu")
FULL.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Pheylalanine',  "Location"), "estimate" = c(coef(modelFULL)$AA[2,1],coef(modelFULL)$AA[2,2], coef(modelFULL)$AA[2,5]) ,
                "std.error" = c(ranef(modelFULL)$AA[2,1], ranef(modelFULL)$AA[2,3], ranef(modelFULL)$AA[2,2]), "group" = 'random', model="Asp")
FULL.ASP <- as_tibble(x)
x<-data.frame("term" = c('Intercept', 'Pheylalanine', "Location"), "estimate" = c(coef(modelFULL)$AA[3,1],coef(modelFULL)$AA[3,2],coef(modelFULL)$AA[3,5]) ,
              "std.error" = c(ranef(modelFULL)$AA[3,1], ranef(modelFULL)$AA[3,3],  ranef(modelFULL)$AA[3,2]), "group" = 'random', model="Pro")
FULL.PRO <- as_tibble(x)
x<-data.frame("term" = c('Intercept', 'Pheylalanine',"Location"), "estimate" = c(coef(modelFULL)$AA[4,1],coef(modelFULL)$AA[4,2],coef(modelFULL)$AA[4,5]) ,
              "std.error" = c(ranef(modelFULL)$AA[4,1], ranef(modelFULL)$AA[4,3], ranef(modelFULL)$AA[4,2]), "group" = 'random', model="Val")
FULL.VAL <- as_tibble(x)
x<-data.frame("term" = c('Intercept', 'Pheylalanine', "Location"), "estimate" = c(coef(modelFULL)$AA[5,1],coef(modelFULL)$AA[5,2], coef(modelFULL)$AA[5,5]) ,
              "std.error" = c(ranef(modelFULL)$AA[5,1], ranef(modelFULL)$AA[5,3],  ranef(modelFULL)$AA[1,2]), "group" = 'random', model="Ala")
FULL.ALA <- as_tibble(x)
FULL.mods <- rbind(FULL.fixed, FULL.GLU, FULL.ASP, FULL.PRO, FULL.VAL, FULL.ALA)



color<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), each=3)
color2<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), times=3)
color3<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), each=4)
color4<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), times=4)
color5<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), each=5)
color6<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), times=5)



FULL.plot <-small_multiple(FULL.mods) +
  theme_bw()+
  geom_point (colour=color5)+
  geom_errorbar (colour=color6, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Full") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1)) 


ENV.plot <-small_multiple(ENV.mods) +
  theme_bw()+
  geom_point (colour=color)+
  geom_errorbar (colour=color2, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Environmental") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1)) 

PREY.plot <-small_multiple(PREY.mods) +
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

NUTRIENT.plot <-small_multiple(NUTRIENT.mods) +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("Coefficient Estimate") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Primary Productivity") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1))


pdf(file="Results/Presentation Figures/NutrientDRAFT.pdf", width=5, height=5)
NUTRIENT.plot
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

pdf(file="Results/Figures/CoefPlot2.pdf", width=11, height=11)
ggarrange(NUTRIENT.plot, PREY.plot, ENV.plot, FULL.plot,rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 2, align= 'hv',
          heights =c(5,5))

dev.off()

###############Residual Plots ##############
library(HLMdiag)
HLMresid(modelFULL, level=1)


pdf(file="Results/Figures/ResidualsDiagnostics.pdf", width=16, height=12)
par(mfrow=c(3,4))

y<- na.omit(cbind(subset(data3, Year>1975 & Year<2012)$TP, subset(data3, Year>1975& Year<2012)$Year))
plot(predict(modelFULL), y[,1], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main=expression(paste("WA ",delta^15, "N Phenylalanine")))
abline(0, 1)

plot(predict(phe.WA.Best), residuals(phe.WA.Best), ylab = "Residuals", xlab = "Predicted",
     main=expression(paste("WA ",delta^15, "N Phenylalanine")),  pch=19, cex=0.5)
abline(0, 0)

y<- na.omit(cbind(subset(data.wa, Year>1950 & Year<2010)$d13C.s, subset(data.wa, Year>1950& Year<2010)$Year))
plot(predict(d13C.WA.Best), y[,1], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main=expression(paste("WA ",delta^13, "C")))
abline(0, 1)
plot(predict(d13C.WA.Best), residuals(d13C.WA.Best), ylab = "Residuals", xlab = "Predicted",
     main=expression(paste("WA ",delta^13, "C")),  pch=19, cex=0.5)
abline(0, 0)




y<- na.omit(cbind(subset(data.goa, Year>1950 & Year<2010)$PHE.mean, subset(data.goa, Year>1950& Year<2010)$Year))
plot(predict(phe.GOA.Best), y[,1], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main=expression(paste("GOA ",delta^15, "N Phenylalanine")))
abline(0, 1)
plot(predict(phe.GOA.Best), residuals(phe.GOA.Best), ylab = "Residuals", xlab = "Predicted",
     main=expression(paste("GOA ",delta^15, "N Phenylalanine")),  pch=19, cex=0.5)
abline(0, 0)

y<- na.omit(cbind(subset(data.goa, Year>1950 & Year<2010)$d13C.s, subset(data.goa, Year>1950& Year<2010)$Year))
plot(predict(d13C.GOA.Best), y[,1], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main=expression(paste("GOA ",delta^13, "C")))
abline(0, 1)
plot(predict(d13C.GOA.Best), residuals(d13C.GOA.Best), ylab = "Residuals", xlab = "Predicted",
     main=expression(paste("GOA ",delta^13, "C")),  pch=19, cex=0.5)
abline(0, 0)




y<- na.omit(cbind(subset(data.ebs, Year>1950 & Year<2010)$PHE.mean, subset(data.ebs, Year>1950& Year<2010)$Year))
plot(predict(phe.EBS.Best), y[,1], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main=expression(paste("EBS ",delta^15, "N Phenylalanine")))
abline(0, 1)
plot(predict(phe.EBS.Best), residuals(phe.EBS.Best), ylab = "Residuals", xlab = "Predicted",
     main=expression(paste("EBS ",delta^15, "N Phenylalanine")),  pch=19, cex=0.5)
abline(0, 0)

y<- na.omit(cbind(subset(data.ebs, Year>1950 & Year<2010)$d13C.s, subset(data.ebs, Year>1950& Year<2010)$Year))
plot(predict(d13C.EBS.Best), y[,1], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main=expression(paste("EBS ",delta^13, "C")))
abline(0, 1)
plot(predict(d13C.EBS.Best), residuals(d13C.EBS.Best), ylab = "Residuals", xlab = "Predicted",
     main=expression(paste("EBS ",delta^13, "C")),  pch=19, cex=0.5)
abline(0, 0)

dev.off()

