library(lme4)
library(AICcmodavg)
data <-  read.csv("Data/Compiled/HierarchicalFiles/HierarchicalData_beta2.csv")


########################    Hier Clim Models       ############################

dataCLIM <- subset(data, AA=="VAL"|AA=="Glu"|AA=="ALA"|AA=="ASP")
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
model.selectionCLIM <- ModelSelection.CLIM(dataCLIM, n, dataCLIM$TP)
model.selectionCLIM
length(model.selectionCLIM[,1])
x<-data.frame(model.selectionCLIM)
subset(x, delAICc<=1.97)
modelENV<- lmer(TP~Location.2+PDO*Col.Dis.high+(1|AA), data=dataCLIM)
summary(modelENV)







########################    Hier PREY Models       ############################

data2 <- subset(data, AA=="Glu"|AA=="ALA"|AA=="Val")
dataPrey <-data2 %>% select(allSmolt, 
                            HakeBiomass,
                            Herring.Biomass,
                            AA,
                            Chinook, 
                            HarborSeal,
                            Coho,
                            Chum,
                           
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
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+(1|AA), data=dataframe)),#2
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
subset(x, delAICc<=5)
PREY<-lmer(TP~Location.2+HakeBiomass+(1|AA), data=dataPrey)
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
    AICc(lmer(y~PDO+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high+Location.2+(1|AA), data=dataframe)),
    
    AICc(lmer(y~PDO+Col.Dis.high+Location.2+(1|AA), data=dataframe)),

    AICc(lmer(y~allSmolt+Col.Dis.high+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~allSmolt*Col.Dis.high+Location.2+(1|AA), data=dataframe)),
    
    AICc(lmer(y~PDO+allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO*allSmolt+Location.2+(1|AA), data=dataframe)),
    
    AICc(lmer(y~PDO+Col.Dis.high+allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~PDO+Col.Dis.high*allSmolt+Location.2+(1|AA), data=dataframe)),
    AICc(lmer(y~Col.Dis.high+PDO*allSmolt+Location.2+(1|AA), data=dataframe))

  )
  
  names <- seq(1,n,1)
  model.names <- c( "Location","PDO", "Smolt","Discharge", "PDO, Dishcarge",
                    "Dis Smolt", "Dis*Smolt", "PDO, Smolt", "PDO*Smolt", 
                    "PDO, dis, smolt", "PDO, dis*smolt", "PDO* smolt, dis")
  
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

modelFULL<- lmer(TP~Location.2+allSmolt*Col.Dis.high+(1|AA), data=dataFull)
FULL<- lmer(TP~Location.2+Col.Dis.high*allSmolt+(1|AA), data=dataFull)

summary(modelFULL)
vif(Full)

##########INTERACTION PLOTS#######

new.DATA.low <- data.frame(
  TP=dataFull$TP,
  Col.Dis.high = dataFull$Col.Dis.high,
  allSmolt = rep(min(dataFull$allSmolt), 300), allSmolt=dataFull$allSmolt,
  Location.2=dataFull$Location.2, Col.Dis.high=dataFull$Col.Dis.high, AA=dataFull$AA
)


new.DATA.high <- data.frame(
  TP=dataFull$TP,
  Col.Dis.high = dataFull$Col.Dis.high,
  allSmolt = rep(max(dataFull$allSmolt), 300), allSmolt=dataFull$allSmolt,
  Location.2=dataFull$Location.2, Col.Dis.high=dataFull$Col.Dis.high, AA=dataFull$AA
)


palette(c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'))
pred<- predict(modelFULL, new.DATA.low, interval = "confidence")
col<-c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1')

pred.low <-data.frame(predict(modelFULL, new.DATA.low, interval = "confidence"))
pred.low.mod<- data.frame(TP=c(pred.low[,1]))

pred.high <-data.frame(predict(modelFULL, new.DATA.high))
pred.high.mod<- data.frame(TP=pred.high[,1])

plot(new.DATA.high$Col.Dis.high, pred.high.mod$TP, col=as.factor(new.DATA.high$AA), pch=16)
legend('topright', legend = levels(pred.high.mod$AA), col = 1:4, cex = 0.8, pch = 1)
plot(new.DATA.low$Col.Dis.high, pred.low.mod$TP,  col=as.factor(new.DATA.low$AA), pch=16)
legend.size<-12

Int.Full <- sjPlot::plot_model(modelFULL, type = "int", col=c(col[1], col[2]))+
  theme_bw() + xlab("Smolt Biomass") + ylab("Predicted Trophic Position") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web", subtitle="Interactions")+
  labs(fill="Smolt Biomass")+
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



Int.Full2 <- sjPlot::plot_model(FULL, type = "int", col=c(col[1], col[2]))+
  theme_bw() + xlab("Discharge") + ylab("Predicted Trophic Position") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Food Web", subtitle="Interactions")+
  labs(fill="Smolt Biomass")+
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
color5<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=3)
color6<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=3)
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


modelPREY <- lmer(TP~Location.2+allSmolt+(1|AA), data=dataPrey)
PREY <- tidy(modelPREY)
x <- data.frame("term" = c('Intercept', 'Location', 'Smolts'), "estimate" = fixef(modelPREY), 
                "std.error" = c(PREY[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
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
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
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


color9<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=5)
color10<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=5)

FULL <- tidy(modelFULL)
x <- data.frame("term" = c('Intercept', 'Location','Smolts', 'Discharge','Discharge:Smolt'), "estimate" = fixef(modelFULL), 
                "std.error" = c(FULL[1:5,3]), "group" = c(rep('fixed', 5)), model="fixed")
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

pdf(file="Results/Figures/HCoefPlot.pdf", width=11, height=12)
ggarrange(ENV.plot, PREY.plot, FULL.plot, Int.Full2,rremove("x.text"), 
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



############### Time Series #################


dataTS<-data %>% select(
  TP,
  AA,
  Location.2,
  Year)




dataTS.coastal<- subset(dataTS, Location.2=="Coastal")
dataTS.coastal<- subset(dataTS.coastal, AA=="VAL"|AA=="Glu"|AA=="ALA"|AA=="PRO")

dataTS.coastal <- dataTS.coastal[complete.cases(dataTS.coastal), ]
dataTS.coastal$AA<- factor(dataTS.coastal$AA, levels = c("Glu","ALA",
                                                             "PRO","VAL"),
                             labels = c("Glutamic Acid", "Alanine", "Proline", "Valine")
)


col<-c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7')
TS.coastal <- ggplot(dataTS.coastal, aes(x = Year, y = TP)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  labs(y="Trophic Position")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  geom_smooth(method="gam",aes(color = AA, alpha=0.5))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(dataTS.coastal$AA))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)
TS.coastal




dataTS.ss<- subset(dataTS, Location.2=="Inland")
dataTS.ss<- subset(dataTS.ss, AA=="VAL"|AA=="Glu"|AA=="ALA"|AA=="PRO")

dataTS.ss <- dataTS.ss[complete.cases(dataTS.ss), ]
dataTS.ss$AA<- factor(dataTS.ss$AA, levels = c("Glu","ALA",
                                                         "PRO","VAL"),
                           labels = c("Glutamic Acid", "Alanine", "Proline", "Valine")
)



TS.ss <- ggplot(dataTS.ss, aes(x = Year, y = TP)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  labs(y="Trophic Position")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  geom_smooth(method="gam",aes(color = AA, alpha=0.5))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(dataTS.ss$AA))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)
TS.ss








resid<-cbind(dataFull,resid=resid(modelFULL))
dataresid<-resid %>% select(
  TP,
  AA,
  Location.2,
  Year,
  resid)




dataresid.coastal<- subset(dataresid, Location.2=="Coastal")
dataresid.coastal<- subset(dataresid.coastal, AA=="VAL"|AA=="Glu"|AA=="ALA"|AA=="PRO")

dataresid.coastal <- dataresid.coastal[complete.cases(dataresid.coastal), ]
dataresid.coastal$AA<- factor(dataresid.coastal$AA, levels = c("Glu","ALA",
                                                         "PRO","VAL"),
                           labels = c("Glutamic Acid", "Alanine", "Proline", "Valine")
)


resid.coastal <- ggplot(dataresid.coastal, aes(x = Year, y = resid)) + 
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
dataresid.ss<- subset(dataresid.ss, AA=="VAL"|AA=="Glu"|AA=="ALA"|AA=="PRO")

dataresid.ss <- dataresid.ss[complete.cases(dataresid.ss), ]
dataresid.ss$AA<- factor(dataresid.ss$AA, levels = c("Glu","ALA",
                                                               "PRO","VAL"),
                              labels = c("Glutamic Acid", "Alanine", "Proline", "Valine")
)


resid.ss <- ggplot(dataresid.ss, aes(x = Year, y = resid)) + 
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

TS.coastal
resid.coastal

TS.ss
resid.ss

############### Location By Sex ###########
males <- subset(data, Sex=='M' &TP>2)
length(males$Sex) #92
females <- subset(data, Sex=='F'&TP>2)
length(females$Sex) #105

Sex <-  subset(data, Sex=='M'|Sex=='F')
summary(lm(TP~Sex, data=Sex))

Sex.c<-subset(Sex, Location.2=="Coastal" & AA!="ASP")
Sex.ss<-subset(Sex, Location.2=="Inland"& AA!="ASP")
sex.GLU.c <- subset(Sex.c, AA=="Glu")
summary(lm(TP~Sex, data=sex.GLU.c))

mean(na.omit(subset(Sex, AA=="VAL")$TP)) #3.9
sd(na.omit(subset(Sex, AA=="VAL")$TP)) #0.77
max(na.omit(subset(Sex, AA=="VAL")$TP)) #3.45
min(na.omit(subset(Sex, AA=="VAL")$TP)) #1.66
vallow <- subset(Sex, AA=="VAL")$TP<3.5
length(vallow[vallow== TRUE])#50
valhigh <- subset(Sex, AA=="VAL")$TP>5
length(valhigh[valhigh== TRUE])#12
length(na.omit(subset(Sex, AA=="VAL")$TP))#108

mean(na.omit(subset(Sex, AA=="Glu")$TP)) #4.5
sd(na.omit(subset(Sex, AA=="Glu")$TP)) #0.44
max(na.omit(subset(Sex, AA=="Glu")$TP)) #5.5
min(na.omit(subset(Sex, AA=="Glu")$TP)) #3.4
glulow <- subset(Sex, AA=="Glu")$TP<3.5
length(glulow[glulow== TRUE])#8
gluhigh <- subset(Sex, AA=="Glu")$TP>5
length(gluhigh[gluhigh== TRUE])#18
length(na.omit(subset(Sex, AA=="Glu")$TP))#109

mean(na.omit(subset(Sex, AA=="ALA")$TP)) #4.1
sd(na.omit(subset(Sex, AA=="ALA")$TP)) #0.77
max(na.omit(subset(Sex, AA=="ALA")$TP)) #5.6
min(na.omit(subset(Sex, AA=="ALA")$TP)) #3.1
alalow <- subset(Sex, AA=="ALA")$TP<3.5
length(alalow[alalow== TRUE])#10
alahigh <- subset(Sex, AA=="ALA")$TP>5
length(alahigh[alahigh== TRUE])#21
length(na.omit(subset(Sex, AA=="ALA")$TP))#108

mean(na.omit(subset(Sex, AA=="PRO")$TP)) #4.6
sd(na.omit(subset(Sex, AA=="PRO")$TP)) #0.81
max(na.omit(subset(Sex, AA=="PRO")$TP)) #6.6
min(na.omit(subset(Sex, AA=="PRO")$TP)) #3.0
prolow <- subset(Sex, AA=="PRO")$TP<3.5
length(prolow[prolow== TRUE])#15
prohigh <- subset(Sex, AA=="PRO")$TP>5
length(prohigh[prohigh== TRUE])#37
length(na.omit(subset(Sex, AA=="PRO")$TP))#107

sex.GLU.ss <- subset(Sex.ss, AA=="Glu")
summary(lm(TP~Sex, data=sex.GLU.ss))
mean(na.omit(sex.GLU.ss$TP))
sd(na.omit(sex.GLU.ss$TP))


col.f <-hcl.colors(5, palette = "Set 3", alpha = NULL, rev = FALSE, fixup = TRUE)
col.m <-hcl.colors(5, palette = "Dark 3", alpha = NULL, rev = FALSE, fixup = TRUE)


male.PRO <- subset(males, AA=='PRO'&Location.2=="Coastal")
females.PRO <- subset(females, AA=="PRO"&Location.2=="Coastal")
male.ALA <- subset(males, AA=='ALA'&Location.2=="Coastal")
females.ALA <- subset(females, AA=="ALA"&Location.2=="Coastal")
male.GLU <- subset(males, AA=='Glu'&Location.2=="Coastal")
females.GLU <- subset(females, AA=="Glu"&Location.2=="Coastal")
male.VAL <- subset(males, AA=='VAL'&Location.2=="Coastal")
females.VAL <- subset(females, AA=="VAL"&Location.2=="Coastal")

male.PRO.s <- subset(males, AA=='PRO'&Location.2=="Inland")
females.PRO.s <- subset(females, AA=="PRO"&Location.2=="Inland")
male.ALA.s <- subset(males, AA=='ALA'&Location.2=="Inland")
females.ALA.s <- subset(females, AA=="ALA"&Location.2=="Inland")
male.GLU.s <- subset(males, AA=='Glu'&Location.2=="Inland")
females.GLU.s <- subset(females, AA=="Glu"&Location.2=="Inland")
male.VAL.s <- subset(males, AA=='VAL'&Location.2=="Inland")
females.VAL.s <- subset(females, AA=="VAL"&Location.2=="Inland")


col.f <-hcl.colors(6, palette = "Set 3", alpha = NULL, rev = FALSE, fixup = TRUE)
col.m <-hcl.colors(6, palette = "Dark 3", alpha = NULL, rev = FALSE, fixup = TRUE)

length(male.GLU$TP)#24
length(females.GLU$TP)#36
length(male.GLU.s$TP)#27
length(females.GLU.s$TP)#22

pdf(file="Results/Figures/IsotopesbySex.pdf", width=10, height=8)

par(mfrow=c(2,1), mar=c(3,5,2,2))
boxplot(list(male.GLU$TP, females.GLU$TP, male.ALA$TP, females.ALA$TP,male.PRO$TP, females.PRO$TP,  male.VAL$TP, females.VAL$TP),
        data=data, at = c(1,2, 4,5, 7,8, 10,11),
        names = c("M", "F", "M", "F","M", "F","M", "F"),
        col=(c(col.m[2], col.f[2], col.m[3], col.f[3], col.m[4], col.f[4],  col.m[5], col.f[5])),
        main="A. Coastal", pch=16, ylim=c(1, 7), ylab= "Trophic Position")
text(1.5,1.5, labels = "Glutamic")
text(1.5,1, labels = "Acid")
text(4.5,1.5, labels = "Alanine")
text(7.5,1.5, labels = "Proline")
text(10.5,1.5, labels = "Valine")

boxplot(list(male.GLU.s$TP, females.GLU.s$TP, male.ALA.s$TP, females.ALA.s$TP,male.PRO.s$TP, females.PRO.s$TP,  male.VAL.s$TP, females.VAL.s$TP),
        data=data, at = c(1,2, 4,5, 7,8, 10,11),
        names = c("M", "F", "M", "F","M", "F","M", "F"),
        col=(c(col.m[2], col.f[2], col.m[3], col.f[3], col.m[4], col.f[4],  col.m[5], col.f[5])),
        main="B. Salish Sea",pch=16, ylim=c(1, 7), ylab= "Trophic Position")
text(1.5,1.5, labels = "Glutamic")
text(1.5,1, labels = "Acid")
text(4.5,1.5, labels = "Alanine")
text(7.5,1.5, labels = "Proline")
text(10.5,1.5, labels = "Valine")
#text(10.5,4, labels = "*", cex=2


dev.off()




pairwise.t.test(Sex.c$TP, pool.sd=FALSE, Sex.c$AA, p.adj = "bonf")
pairwise.t.test(Sex.ss$TP, pool.sd=FALSE, Sex.ss$AA, p.adj = "bonf")
summary(lm(TP~AA*Sex, data=Sex.c))
summary(lm(TP~AA*Sex, data=Sex.ss))

############### Location By Length ###########
Length <-data %>% select(TP,
                              AA,
                               Length,
                               Location.2)

Length <- Length[complete.cases(Length), ]



Length.c <- Length%>% filter(Location.2 == "Coastal"&AA!="ASP")
Length.c$AA<- factor(Length.c$AA, levels = c("Glu","ALA",
                                                         "PRO","VAL"),
                           labels = c("Glutamic Acid", "Alanine", "Proline", "Valine")
)
summary(lm(TP~AA*Length, data=Length.c))
fit.c<-lm(TP~AA*Length, data=Length.c)

Coastal.Length <- qplot(Length, TP, data=Length.c, colour=AA)+
  ggtitle("A. Coastal")+
  theme_bw()+
  labs(y="Trophic Position", x="Standard Length (cm)")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  #geom_smooth(method="lm", aes(color = AA, alpha=0.5))+
  geom_line(data = fortify(fit.c), aes(x = Length, y = .fitted))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(Length.c$AA))+
  theme(plot.title = element_text(hjust = 0.5),strip.background =element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"
        ))+
  guides(colour=FALSE, alpha=FALSE)
Coastal.Length 


Length.ss <- Length%>% filter(Location.2 == "Inland"&AA!="ASP")
Length.ss$AA<- factor(Length.ss$AA, levels = c("Glu","ALA",
                                             "PRO","VAL"),
                     labels = c("Glutamic Acid", "Alanine", "Proline", "Valine")
)
summary(lm(TP~AA*Length, data=Length.ss))
fit.ss<-lm(TP~AA*Length, data=Length.ss)

SalishSea.Length <- qplot(Length, TP, data=Length.ss, colour=AA)+
  ggtitle("B. Salish Sea")+
  theme_bw()+
  labs(y="", x="Standard Length (cm)")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  #geom_smooth(method="lm", aes(color = AA, alpha=0.5))+
  geom_line(data = fortify(fit.ss), aes(x = Length, y = .fitted))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(Length.ss$AA))+
  theme(plot.title = element_text(hjust = 0.5), strip.background =element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"
        ))+
  guides(colour=FALSE, alpha=FALSE)
SalishSea.Length 

pdf(file="Results/Figures/Lengthplot.pdf", width=8, height=4)
ggarrange(Coastal.Length, SalishSea.Length, rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 1, align= 'hv')

dev.off()

