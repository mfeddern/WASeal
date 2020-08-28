library(lme4)
library(AICcmodavg)
data <-  read.csv("Data/Compiled/HierarchicalData2.csv")


########################    Hier Clim Models       ############################

dataCLIM <- subset(data, AA=="VAL"|AA=="Glu"|AA=="PRO"|AA=="ALA")
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
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe)),
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe)), 
                      AICc(lmer(y~PDO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~NPGO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~MEI+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~PDO+NPGO+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~PDO+MEI+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~MEI+NPGO+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      AICc(lmer(y~MEI+NPGO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#7
                      AICc(lmer(y~PDO+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#8
                      AICc(lmer(y~PDO+MEI+NPGO+Location.2+(1|AA), data=dataframe)),#9
                      AICc(lmer(y~PDO+NPGO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      
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
                      
                      AICc(lmer(y~Location.2+(1|AA)+Col.Dis.high, data=dataframe)), 
                      AICc(lmer(y~Col.Dis.high+PDO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+NPGO+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+MEI+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+PDO+NPGO+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Col.Dis.high+PDO+MEI+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~Col.Dis.high+MEI+NPGO+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe)),#6
                      
                      AICc(lmer(y~Col.Dis.high+WA.SST.Su+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe)),#6
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2+(1|AA), data=dataframe))#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Null","Location", "PDO", "NPGO", "MEI", "Upwelling (Sp)","1. PDO, NPGO", "2. PDO, Upwelling (Sp)", "3. PDO, MEI", "4. NPGO, Upwelling (Sp)",
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
model.selectionCLIM <- ModelSelection.CLIM(dataCLIM, n, dataCLIM$TP.norm)
model.selectionCLIM
length(model.selectionCLIM[,1])
x<-data.frame(model.selectionCLIM)
subset(x, delAICc<=1.97)
modelENV<- lmer(TP~Location.2+Col.Dis.high+PDO+(1|AA), data=dataCLIM)
summary(modelENV)







########################    Hier PREY Models       ############################

data2 <- subset(data, AA=="Glu"|AA=="PRO"|AA=="ALA"|AA=="VAL")
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


model.selection.PREY <- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe)), 
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Chinook+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~allSmolt+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),#7
                      AICc(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#8
                      AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#9
                      AICc(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      
                      AICc(lmer(y~Chum+Location.2+(1|AA), data=dataframe)), 
                      AICc(lmer(y~Chum+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Chum+Chinook+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Chum+allSmolt+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Chum+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Chum+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~Chum+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Chum+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~Chum+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~Chum+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~Chum+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      
                      AICc(lmer(y~Coho+Location.2+(1|AA), data=dataframe)), 
                      AICc(lmer(y~Coho+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Coho+Chinook+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Coho+allSmolt+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Coho+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~Coho+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~Coho+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Coho+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~Coho+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~Coho+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~Coho+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      
                      AICc(lmer(y~HarborSeal+Location.2+(1|AA), data=dataframe)), 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Chinook+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+allSmolt+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe)),#3
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+allSmolt+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      
                      AICc(lmer(y~HarborSeal+Chum+Location.2+(1|AA), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+Coho+Location.2+(1|AA), data=dataframe)),#5
                      AICc(lmer(y~Chum+Coho+Location.2+(1|AA), data=dataframe)),#6
                      AICc(lmer(y~Chum+Coho+HarborSeal+Location.2+(1|AA), data=dataframe)),#6
                      
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+allSmolt+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+Coho+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+HarborSeal+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+Chum+Location.2+(1|AA), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+allSmolt+HarborSeal+Location.2+(1|AA), data=dataframe))#2
                      
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Null","Location", "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
                   "5. Chinook, Hatch", "6. Hatch, Hake", "7. Chinook, Hake, Hatch", "8. Herring, Hatch, Hake", "9. Herring, Hatch, Chinook", 
                   "11. Herring Hake Chinook",
                   
                   "Chum, Location", "Chum, Herring", "Chum, Chinook", "Chum, Hatch", "Chum, Hake","1. Chum,  Herring, Chinook", "2. Chum, Herring, Hake", 
                   "3. Chum, Herring, Hatch", "4. Chum, Chinook, Hake","5. Chum, Chinook, Hatch", "6. Chum, Hatch, Hake",
                   
                   "Coho, Location", "Coho, Herring", "Coho, Chinook", "Coho, Hatch", "Coho, Hake","1. Coho, Herring, Chinook", "2.Coho,  Herring, Hake",
                   "3. Coho, Herring, Hatch", "4.Coho,  Chinook, Hake","5.Coho,  Chinook, Hatch", "6. Coho, Hatch, Hake",
                   
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook", "Harbor Seal, Hatch", "Hake","1. Harbor Seal, Herring, Chinook",
                   "2.Harbor Seal,  Herring, Hake", "3.Harbor Seal,  Herring, Hatch", "4.Harbor Seal,  Chinook, Hake","5. Harbor Seal, Chinook, Hatch", 
                   "6. Harbor Seal, Hatch, Hake",
                   
                   "Harbor Seal, Chum", " Harbor Seal Coho", "Chum, Coho", "Chum, Coho, Harbor Seal",
                   "Int", "Int smolt", "Int harbor seal", "Int Coho", "Int Chum", "Int smolt seal")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}


model.selectionPREY <- model.selection.PREY(dataPrey, n, dataPrey$TP)
length(model.selectionPREY[,1])
model.selectionPREY
x<-data.frame(model.selectionPREY)
subset(x, delAICc<=1.97)
PREY<-lmer(TP~Location.2+HakeBiomass*Herring.Biomass+allSmolt+(1|AA), data=dataPrey)
summary(PREY)
vif(PREY)


########################     Hierarchical Nutrient Models       ############################
dataNut <- subset(data,  AA=="Glu"|AA=="PRO"|AA=="ALA"|AA=="VAL")
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
    AICc(lmer(y~(1|AA), data=dataframe)),
    AICc(lmer(y~Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PHE.norm+d13C.norm+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PHE.norm+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~d13C.norm+Location.2+(1|AA), data=dataframe))
  )
  
  names <- seq(1,n,1)
  model.names <- c( "Null","Location","Phe, 13C", "PHE", "13C")
  
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
model.selectionNUTRIENT <- ModelSelection.WANut(dataNut, n, dataNut$TP)
model.selectionNUTRIENT
x<-data.frame(model.selectionNUTRIENT)
subset(x, delAICc<=1.97)
summary(lmer(TP~PHE.norm+Location.2+(1|AA), data=dataNut))
modelNUTRIENT<-lmer(TP~PHE.norm+Location.2+(1|AA), data=dataNut)
########################     Hierarchical Full Models       ############################
dataFull <- subset(data,  AA=="Glu"|AA=="PRO"|AA=="ALA"|AA=="VAL")
dataFull <-dataFull %>% select(PHE.mean,
                             PHE.norm,
                             AA,
                             allSmolt, 
                             HakeBiomass,
                             Herring.Biomass,
                             d13C.norm,
                             TP.norm,
                             TP,
                             Col.Dis.high,
                             PDO,
                             Year,
                             Sample.ID,
                             Location.2)
dataFull <- dataFull[complete.cases(dataFull), ]



ModelSelection.WAFull <- function(dataframe,n, y) {
  
  aic.output <- rbind(
    AICc(lmer(y~Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO+Col.Dis.high+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~allSmolt+Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO+Location.2+(1|AA), data=dataframe)),
    
    AICc(lmer(y~PDO+allSmolt+Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO+Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO*allSmolt+Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
    
    
    AICc(lmer(y~PDO+allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO*allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO*Herring.Biomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO+HakeBiomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO*HakeBiomass+Location.2+(1|AA), data=dataframe)),

    AICc(lmer(y~Col.Dis.high+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high+allSmolt+Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high+Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high*allSmolt+Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe)),
    
    
    AICc(lmer(y~Col.Dis.high+allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high*allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high+Herring.Biomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high*Herring.Biomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high+HakeBiomass+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high*HakeBiomass+Location.2+(1|AA), data=dataframe))
    
  )
  
  names <- seq(1,n,1)
  model.names <- c( "Location","Environment", "Prey","PDO", "PDO, Prey","PDO, Herring*Hake",
                    "PDO*Smolt Hake*Herring", "PDO, Smolt", "PDO*Smolt", "PDO, Herring", 
                    "PDO*Herring", "PDO, Hake", "PDO*Hake",
                    
                    
                    "Col", "Col, Prey","Col, Herring*Hake","Col*Smolt Hake*Herring", "Col, Smolt",
                    "Col*Smolt", "Col, Herring", 
                    "Col*Herring", "Col, Hake", "Col*Hake")
  
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
model.selectionFULL <- ModelSelection.WAFull(dataFull, n, dataFull$TP)
model.selectionFULL
x<-data.frame(model.selectionFULL)
subset(x, delAICc<=1.97)

modelFULL<- lmer(TP~Location.2+allSmolt+HakeBiomass*Herring.Biomass+Col.Dis.high+(1|AA), data=dataFull)

summary(Full)
vif(Full)

##########INTERACTION PLOTS#######

new.DATA.low <- data.frame(
  TP=dataFull$TP,
  HakeBiomass = dataFull$HakeBiomass,
  Herring.Biomass = rep(-1, 300), allSmolt=dataFull$allSmolt,
  Location.2=dataFull$Location.2, Col.Dis.high=dataFull$Col.Dis.high, AA=dataFull$AA
)


new.DATA.high <- data.frame(
  TP=dataFull$TP,
  HakeBiomass = dataFull$HakeBiomass,
  Herring.Biomass = rep(2.2, 300), allSmolt=dataFull$allSmolt,
  Location.2=dataFull$Location.2, Col.Dis.high=dataFull$Col.Dis.high, AA=dataFull$AA
)


palette(c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'))
pred<- predict(PREY, new.DATA.low, interval = "confidence")
col<-c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1')

pred.low <-data.frame(predict(PREY, new.DATA.low, interval = "confidence"))
pred.low.mod<- data.frame(TP=c(pred.low[,1]))

pred.high <-data.frame(predict(PREY, new.DATA.high))
pred.high.mod<- data.frame(TP=pred.high[,1])

plot(new.DATA.high$HakeBiomass, pred.high.mod$TP, col=as.factor(new.DATA.high$AA), pch=16)
legend('topright', legend = levels(pred.high.mod$AA), col = 1:4, cex = 0.8, pch = 1)
plot(new.DATA.low$HakeBiomass, pred.low.mod$TP,  col=as.factor(new.DATA.low$AA), pch=16)
legend.size<-12

Int.Prey <- sjPlot::plot_model(PREY, type = "int", col=c(col[1], col[2]))+
  theme_bw() + xlab("Hake Biomass") + ylab("Predicted Trophic Position") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web", subtitle="Interactions")+
  labs(fill="Herring Biomass")+
 # scale_color_manual(values= c("#CCA65A", "#7EBA68"),name="Hake Biomass", labels = c("Low", "High"))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title=element_text(size=13),
        legend.title.align = 0,
        legend.title=element_text(size=legend.size),
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))



Int.Full <- sjPlot::plot_model(modelFULL, type = "int", col=c(col[1], col[2]))+
  theme_bw() + xlab("Hake Biomass") + ylab("Predicted Trophic Position") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web", subtitle="Interactions")+
  labs(fill="Herring Biomass")+
  # scale_color_manual(values= c("#CCA65A", "#7EBA68"),name="Hake Biomass", labels = c("Low", "High"))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0.5, 0.25, 0.25, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title=element_text(size=13),
        legend.title.align = 0,
        legend.title=element_text(size=legend.size),
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))



new.DATA.low <- data.frame(
  TP=dataFull$TP,
  HakeBiomass = rep(-1.5, 300),
  Herring.Biomass = dataFull$Herring.Biomass, allSmolt=dataFull$allSmolt,
  Location.2=dataFull$Location.2, Col.Dis.high=dataFull$Col.Dis.high, AA=dataFull$AA
)


new.DATA.high <- data.frame(
  TP=dataFull$TP,
  HakeBiomass = rep(1.5, 300),
  Herring.Biomass = dataFull$Herring.Biomass, allSmolt=dataFull$allSmolt,
  Location.2=dataFull$Location.2, Col.Dis.high=dataFull$Col.Dis.high, AA=dataFull$AA
)

pred<- predict(Full, new.DATA.low, interval = "confidence")
pred.low <-data.frame(predict(Full, new.DATA.low, interval = "confidence"))
pred.low.mod<- data.frame(TP=c(pred.low[,1]))

pred.high <-data.frame(predict(Full, new.DATA.high))
pred.high.mod<- data.frame(TP=pred.high[,1])

plot(new.DATA.high$Herring.Biomass, pred.high.mod$TP, col=as.factor(new.DATA.high$AA), pch=16)
plot(new.DATA.low$Herring.Biomass, pred.low.mod$TP,  col=as.factor(new.DATA.low$AA), pch=16)

########################     Plots#######################
library(dotwhisker)
library(broom)
library(dplyr)
library(colorspace)
library(ggpubr)

color<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=4)
color2<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=4)
color3<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=3)
color4<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=3)
color5<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=6)
color6<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=6)
color7<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), each=5)
color8<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), times=5)



ENV <- tidy(modelENV)
x <- data.frame("term" = c('Intercept', 'Location', 'Discharge', 'PDO'), "estimate" = fixef(modelENV), "std.error" = c(ENV[1:4,3]), "group" = c(rep('fixed', 4)), model="fixed")
ENV.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV)$AA[1,1]), "std.error" = c(ranef(modelENV)$AA[1,1]), "group" = 'random', model="Ala")
ENV.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV)$AA[2,1]), "std.error" = c(ranef(modelENV)$AA[2,1]), "group" = 'random', model="Glu")
ENV.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV)$AA[3,1]), "std.error" = c(ranef(modelENV)$AA[3,1]), "group" = 'random', model="Pro")
ENV.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelENV)$AA[4,1]), "std.error" = c(ranef(modelENV)$AA[4,1]), "group" = 'random', model="Val")
ENV.VAL <- as_tibble(x)

ENV.mods <- rbind(ENV.fixed, ENV.GLU, ENV.ALA, ENV.PRO, ENV.VAL)

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
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10)) 


NUTRIENT <- tidy(modelNUTRIENT)
x <- data.frame("term" = c('Intercept','Phenylalanine','Location'), "estimate" = fixef(modelNUTRIENT),
                "std.error" = c(NUTRIENT[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
NUTRIENT.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelNUTRIENT)$AA[1,1]) ,
                "std.error" = c(ranef(modelNUTRIENT)$AA[1,1]), "group" = 'random', model="Ala")
NUTRIENT.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelNUTRIENT)$AA[2,1]) ,
                "std.error" = c(ranef(modelNUTRIENT)$AA[2,1]), "group" = 'random', model="Glu")
NUTRIENT.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelNUTRIENT)$AA[3,1]) ,
                "std.error" = c(ranef(modelNUTRIENT)$AA[3,1]), "group" = 'random', model="Pro")
NUTRIENT.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelNUTRIENT)$AA[4,1]) ,
                "std.error" = c(ranef(modelNUTRIENT)$AA[4,1]), "group" = 'random', model="Val")
NUTRIENT.VAL <- as_tibble(x)
NUTRIENT.mods <- rbind(NUTRIENT.fixed, NUTRIENT.GLU,NUTRIENT.ALA,  NUTRIENT.PRO, NUTRIENT.VAL)

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
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))


modelPREY <- lmer(TP~Location.2+Herring.Biomass*HakeBiomass+allSmolt+(1|AA), data=dataPrey)
PREY <- tidy(modelPREY)
x <- data.frame("term" = c('Intercept', 'Location',  'Herring', 'Hake', 'Smolts', 'Hake:Herring'), "estimate" = fixef(modelPREY), 
                "std.error" = c(PREY[1:6,3]), "group" = c(rep('fixed', 6)), model="fixed")
PREY.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[1,1]) ,
                "std.error" = c(ranef(modelPREY)$AA[1,1]), "group" = 'random', model="Ala")
PREY.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[2,1]),
                "std.error" = c(ranef(modelPREY)$AA[2,1]), "group" = 'random', model="Glu")
PREY.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[3,1]), 
                "std.error" = c(ranef(modelPREY)$AA[3,1]), "group" = 'random', model="Pro")
PREY.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelPREY)$AA[4,1]), 
                "std.error" = c(ranef(modelPREY)$AA[4,1]), "group" = 'random', model="Val")
PREY.VAL <- as_tibble(x)
PREY.mods <- rbind(PREY.fixed, PREY.GLU, PREY.ALA,  PREY.PRO, PREY.VAL)


PREY.plot <-small_multiple(PREY.mods) +
  theme_bw()+
  geom_point (colour=color5)+
  geom_errorbar (colour=color6, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web") +
 # scale_y_discrete("", breaks=waiver(), labels=waiver())+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))


color9<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=7)
color10<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=7)

FULL <- tidy(modelFULL)
x <- data.frame("term" = c('Intercept', 'Location','Smolts',  'Herring', 'Hake','Discharge','Hake:Herring'), "estimate" = fixef(modelFULL), 
                "std.error" = c(FULL[1:7,3]), "group" = c(rep('fixed', 7)), model="fixed")
FULL.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelFULL)$AA[1,1]) ,
                "std.error" = c(ranef(modelFULL)$AA[1,1]), "group" = 'random', model="Ala")
FULL.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelFULL)$AA[2,1]),
                "std.error" = c(ranef(modelFULL)$AA[2,1]), "group" = 'random', model="Glu")
FULL.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelFULL)$AA[3,1]), 
                "std.error" = c(ranef(modelFULL)$AA[3,1]), "group" = 'random', model="Pro")
FULL.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept'), "estimate" = c(coef(modelFULL)$AA[4,1]), 
                "std.error" = c(ranef(modelFULL)$AA[4,1]), "group" = 'random', model="Val")
FULL.VAL <- as_tibble(x)
FULL.mods <- rbind(FULL.fixed, FULL.GLU, FULL.ALA,  FULL.PRO, FULL.VAL)




FULL.plot <-small_multiple(FULL.mods) +
  theme_bw()+
  geom_point (colour=color9)+
  geom_errorbar (colour=color10, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Full") +
  ylab("Coefficient Estimate") +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))


color3<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=6)
color4<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=6)




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

pdf(file="Results/Figures/HCoefPlot2.pdf", width=11, height=12)
ggarrange(NUTRIENT.plot, ENV.plot, PREY.plot, Int.Prey,rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 2, align= 'hv',
          heights =c(4,6))

dev.off()

pdf(file="Results/Figures/HCoefPlot2.pdf", width=11, height=12)
ggarrange(NUTRIENT.plot, ENV.plot, PREY.plot, Int.Prey,rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 2, align= 'hv',
          heights =c(4,6))

dev.off()


pdf(file="Results/Figures/HCoefPlotFULL2.pdf", width=11, height=7)
ggarrange(FULL.plot, Int.Full,rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 1, align= 'hv',
          heights =c(6))

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

