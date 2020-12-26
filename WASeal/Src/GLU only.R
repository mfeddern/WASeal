rm(list = ls())

library(lme4)
library(AICcmodavg)

dataALL <-  read.csv("Data/Compiled/HierarchicalData2.csv")
data<-subset(subset(dataALL, Location.2=="Inland"|Location.2=="Coastal"), beta==1& eq==2)
data$TP
subset(dataALL, AA=="GLU")

########################  
########################    Hier Clim2 Models       ############################

dataCLIM <- subset(data, AA=="GLU")
#dataCLIM <- subset(data, Year>=1960&Year<=2008)# & AA=="Glu"|AA=="PRO"|AA=="ALA")

dataCLIM <-dataCLIM %>% select(MEI, 
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



dataCLIM <- dataCLIM[complete.cases(dataCLIM), ]
ModelSelection.CLIM<- function(dataframe,n, y){
  
  aic.output <- rbind( #1
                      AICc(lm(y~Location.2 , data=dataframe)), #2
                      AICc(lm(y~PDO+Location.2 , data=dataframe)),#3
                      AICc(lm(y~NPGO+Location.2 , data=dataframe)),#4
                      AICc(lm(y~MEI+Location.2 , data=dataframe)),#5
                      AICc(lm(y~UpInAn.45.Spring+Location.2 , data=dataframe)),#6
                      AICc(lm(y~PDO+NPGO+Location.2 , data=dataframe)),
                      AICc(lm(y~PDO+UpInAn.45.Spring+Location.2 , data=dataframe)),#7
                      #AICc(lm(y~PDO+MEI+Location.2 , data=dataframe)),
                      AICc(lm(y~UpInAn.45.Spring+NPGO+Location.2 , data=dataframe)),#8
                      #AICc(lm(y~MEI+NPGO+Location.2 , data=dataframe)),#5
                      AICc(lm(y~MEI+UpInAn.45.Spring+Location.2 , data=dataframe)),#9
                      #AICc(lm(y~MEI+NPGO+UpInAn.45.Spring+Location.2 , data=dataframe)),#7
                      #AICc(lm(y~PDO+MEI+UpInAn.45.Spring+Location.2 , data=dataframe)),#8
                      #AICc(lm(y~PDO+MEI+NPGO+Location.2 , data=dataframe)),#9
                      #AICc(lm(y~PDO+NPGO+UpInAn.45.Spring+Location.2 , data=dataframe)),
                      
                      AICc(lm(y~Location.2 +WA.SST.Su, data=dataframe)), #10
                      AICc(lm(y~WA.SST.Su+PDO+Location.2 , data=dataframe)),#11
                      AICc(lm(y~WA.SST.Su+NPGO+Location.2 , data=dataframe)),#12
                      AICc(lm(y~WA.SST.Su+MEI+Location.2 , data=dataframe)),#13
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Spring+Location.2 , data=dataframe)),#14
                      # AIC(lm(y~WA.SST.Su+PDO+NPGO+Location.2 , data=dataframe)),
                      AICc(lm(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2 , data=dataframe)),#15
                      #AIC(lm(y~WA.SST.Su+PDO+MEI+Location.2 , data=dataframe)),
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2 , data=dataframe)),#16
                      #AIC(lm(y~WA.SST.Su+MEI+NPGO+Location.2 , data=dataframe)),
                      AICc(lm(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2 , data=dataframe)),#17
                      
                      AICc(lm(y~Location.2 +UpInAn.45.Summer, data=dataframe)), #18
                      AICc(lm(y~UpInAn.45.Summer+PDO+Location.2 , data=dataframe)),#19
                      AICc(lm(y~UpInAn.45.Summer+NPGO+Location.2 , data=dataframe)),#20
                      AICc(lm(y~UpInAn.45.Summer+MEI+Location.2 , data=dataframe)),#21
                      AICc(lm(y~UpInAn.45.Summer+UpInAn.45.Spring+Location.2 , data=dataframe)),#22
                      # AIC(lm(y~UpInAn.45.Summer+PDO+NPGO+Location.2 , data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2 , data=dataframe)),
                      #AIC(lm(y~UpInAn.45.Summer+PDO+MEI+Location.2 , data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2 , data=dataframe)),#23
                      #AIC(lm(y~UpInAn.45.Summer+MEI+NPGO+Location.2 , data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2 , data=dataframe)),#24
                      
                      AICc(lm(y~Location.2 +Col.Dis.high, data=dataframe)), #25
                      AICc(lm(y~Col.Dis.high+PDO+Location.2 , data=dataframe)),#26
                      AICc(lm(y~Col.Dis.high+NPGO+Location.2 , data=dataframe)),#27
                      AICc(lm(y~Col.Dis.high+MEI+Location.2 , data=dataframe)),#28
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Spring+Location.2 , data=dataframe)),#29
                      # AICc(lm(y~Col.Dis.high+PDO+NPGO+Location.2 , data=dataframe)),
                      AICc(lm(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2 , data=dataframe)),#30
                      #AICc(lm(y~Col.Dis.high+PDO+MEI+Location.2 , data=dataframe)),
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2 , data=dataframe)),#31
                      # AICc(lm(y~Col.Dis.high+MEI+NPGO+Location.2 , data=dataframe)),
                      AICc(lm(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2 , data=dataframe)),#32
                      
                      AICc(lm(y~Col.Dis.high+WA.SST.Su+Location.2 , data=dataframe)),#33
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Summer+Location.2 , data=dataframe)),#34
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Summer+Location.2 , data=dataframe)),#35
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2 , data=dataframe)),
                      AICc(lm(y~MEI*Location.2 , data=dataframe))#36
                      #36
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("2. Location", "3. PDO", "4. NPGO", "5. MEI", "6. Upwelling (Sp)","NPGO,PDO", "7. PDO, Upwelling (Sp)",  
                   "8. NPGO, Upwelling (Sp)","9. MEI, Upwelling (Sp)",  
                   
                   "10. WA.SST.Su, Location", "11. WA.SST.Su, PDO", "12. WA.SST.Su, NPGO", "13. WA.SST.Su, MEI",
                   "14. WA.SST.Su, Upwelling (Sp)", "15. WA.SST.Su, PDO, Upwelling (Sp)", 
                   "16. WA.SST.Su, NPGO, Upwelling (Sp)", "17. WA.SST.Su, MEI, Upwelling (Sp)",
                   
                   "18. UpInAn.45.Summer, Location", "19. UpInAn.45.Summer, PDO", "20. UpInAn.45.Summer, NPGO", 
                   "21. UpInAn.45.Summer, MEI", "22. UpInAn.45.Summer, Upwelling (Sp)",
                   "23.UpInAn.45.Summer,  PDO, Upwelling (Sp)",
                   "24. UpInAn.45.Summer,  NPGO, Upwelling (Sp)", "25. UpInAn.45.Summer, MEI, Upwelling (Sp)",
                   
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
n<- 17
model.selectionCLIM <- ModelSelection.CLIM(dataCLIM, n, dataCLIM$TP)
model.selectionCLIM
length(model.selectionCLIM[,1])
x<-data.frame(model.selectionCLIM)
subset(x, delAICc<=1.97)
modelENV<- lm(TP~Location.2+MEI , data=dataCLIM)
summary(modelENV)
VIF(lm(TP~Location.2+NPGO+MEI , data=dataCLIM))
vif(lm(TP~Location.2+NPGO+MEI , data=dataCLIM))


########################    Hier PREY Models       ############################

data2 <- subset(data, AA=="GLU")
dataPrey <-data2 %>% select(allSmolt, 
                            HakeBiomass,
                            Herring.Biomass,
                            AA,
                            Chinook, 
                            HarborSeal,
                            Coho,
                            Chum,
                            MEI,
                            TP,
                            Year,
                            Sample.ID,
                            Location.2,
                            eq,
                            beta)
dataPrey <- dataPrey[complete.cases(dataPrey), ]


model.selection.PREY <- function(dataframe,n, y) {
  
  aic.output <- rbind( 
                      AICc(lm(y~Location.2 , data=dataframe)),
                      AICc(lm(y~Herring.Biomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Chinook+Location.2 , data=dataframe)),
                      AICc(lm(y~allSmolt+Location.2 , data=dataframe)),
                      AICc(lm(y~HakeBiomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Herring.Biomass+Chinook+Location.2 , data=dataframe)),#1
                      AICc(lm(y~Herring.Biomass+HakeBiomass , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass+allSmolt+Location.2 , data=dataframe)),#3
                      AICc(lm(y~HakeBiomass+Chinook+Location.2 , data=dataframe)),#4
                      AICc(lm(y~allSmolt+Chinook+Location.2 , data=dataframe)),#5
                      AICc(lm(y~allSmolt+HakeBiomass+Location.2 , data=dataframe)),#6
                      AICc(lm(y~allSmolt+Chinook+HakeBiomass+Location.2 , data=dataframe)),#7
                      AICc(lm(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2 , data=dataframe)),#8
                      AICc(lm(y~Herring.Biomass+allSmolt+Chinook+Location.2 , data=dataframe)),#9
                      AICc(lm(y~Herring.Biomass+Chinook+HakeBiomass+Location.2 , data=dataframe)),
                      
                      AICc(lm(y~Chum+Location.2 , data=dataframe)), 
                      AICc(lm(y~Chum+Herring.Biomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Chum+Chinook+Location.2 , data=dataframe)),
                      AICc(lm(y~Chum+allSmolt+Location.2 , data=dataframe)),
                      AICc(lm(y~Chum+HakeBiomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Chum+Herring.Biomass+Chinook+Location.2 , data=dataframe)),#1
                      AICc(lm(y~Chum+Herring.Biomass+HakeBiomass+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Chum+Herring.Biomass+allSmolt+Location.2 , data=dataframe)),#3
                      AICc(lm(y~Chum+HakeBiomass+Chinook+Location.2 , data=dataframe)),#4
                      AICc(lm(y~Chum+allSmolt+Chinook+Location.2 , data=dataframe)),#5
                      AICc(lm(y~Chum+allSmolt+HakeBiomass+Location.2 , data=dataframe)),#6
                      
                      AICc(lm(y~Coho+Location.2 , data=dataframe)), 
                      AICc(lm(y~Coho+Herring.Biomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Coho+Chinook+Location.2 , data=dataframe)),
                      AICc(lm(y~Coho+allSmolt+Location.2 , data=dataframe)),
                      AICc(lm(y~Coho+HakeBiomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Coho+Herring.Biomass+Chinook+Location.2 , data=dataframe)),#1
                      AICc(lm(y~Coho+Herring.Biomass+HakeBiomass+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Coho+Herring.Biomass+allSmolt+Location.2 , data=dataframe)),#3
                      AICc(lm(y~Coho+HakeBiomass+Chinook+Location.2 , data=dataframe)),#4
                      AICc(lm(y~Coho+allSmolt+Chinook+Location.2 , data=dataframe)),#5
                      AICc(lm(y~Coho+allSmolt+HakeBiomass+Location.2 , data=dataframe)),#6
                      
                      AICc(lm(y~HarborSeal+Location.2 , data=dataframe)), 
                      AICc(lm(y~HarborSeal+Herring.Biomass+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+Chinook+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+allSmolt+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+HakeBiomass+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+Herring.Biomass+Chinook+Location.2 , data=dataframe)),#1
                      AICc(lm(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2 , data=dataframe)),#2
                      AICc(lm(y~HarborSeal+Herring.Biomass+allSmolt+Location.2 , data=dataframe)),#3
                      AICc(lm(y~HarborSeal+HakeBiomass+Chinook+Location.2 , data=dataframe)),#4
                      AICc(lm(y~HarborSeal+allSmolt+Chinook+Location.2 , data=dataframe)),#5
                      AICc(lm(y~HarborSeal+allSmolt+HakeBiomass+Location.2 , data=dataframe)),#6
                      
                      AICc(lm(y~HarborSeal+Chum+Location.2 , data=dataframe)),#4
                      AICc(lm(y~HarborSeal+Coho+Location.2 , data=dataframe)),#5
                      AICc(lm(y~Chum+Coho+Location.2 , data=dataframe)),#6
                      AICc(lm(y~Chum+Coho+HarborSeal+Location.2 , data=dataframe)),#6
                      
                      AICc(lm(y~allSmolt*HakeBiomass+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*allSmolt+Location.2 , data=dataframe)),#2
                      AICc(lm(y~allSmolt*HarborSeal+Coho+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HakeBiomass+HarborSeal+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HakeBiomass+Chum+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HakeBiomass+allSmolt+HarborSeal+Location.2 , data=dataframe))#2
                      
                      
                      
                      
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


model.selection.PREY2 <- function(dataframe,n, y) {
  
  aic.output <- rbind( 
                      AICc(lm(y~Location.2 , data=dataframe)),
                      AICc(lm(y~Herring.Biomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Chinook+Location.2 , data=dataframe)),
                      AICc(lm(y~allSmolt+Location.2 , data=dataframe)),
                      AICc(lm(y~HakeBiomass+Location.2 , data=dataframe)),
                      AICc(lm(y~Herring.Biomass+Chinook+Location.2 , data=dataframe)),#1
                      AICc(lm(y~Herring.Biomass+HakeBiomass , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass+allSmolt+Location.2 , data=dataframe)),#3
                      AICc(lm(y~HakeBiomass+Chinook+Location.2 , data=dataframe)),#4
                      AICc(lm(y~allSmolt+Chinook+Location.2 , data=dataframe)),#5
                      AICc(lm(y~allSmolt+HakeBiomass+Location.2 , data=dataframe)),#6
                      AICc(lm(y~allSmolt+Chinook+HakeBiomass+Location.2 , data=dataframe)),#7
                      AICc(lm(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2 , data=dataframe)),#8
                      AICc(lm(y~Herring.Biomass+allSmolt+Chinook+Location.2 , data=dataframe)),#9
                      AICc(lm(y~Herring.Biomass+Chinook+HakeBiomass+Location.2 , data=dataframe)),
                      
                      
                      AICc(lm(y~HarborSeal+Location.2 , data=dataframe)), 
                      AICc(lm(y~HarborSeal+Herring.Biomass+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+Chinook+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+allSmolt+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+HakeBiomass+Location.2 , data=dataframe)),
                      AICc(lm(y~HarborSeal+Herring.Biomass+Chinook+Location.2 , data=dataframe)),#1
                      AICc(lm(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2 , data=dataframe)),#2
                      AICc(lm(y~HarborSeal+Herring.Biomass+allSmolt+Location.2 , data=dataframe)),#3
                      AICc(lm(y~HarborSeal+HakeBiomass+Chinook+Location.2 , data=dataframe)),#4
                      #AICc(lm(y~HarborSeal+allSmolt+Chinook+Location.2 , data=dataframe)),#5
                      AICc(lm(y~HarborSeal+allSmolt+HakeBiomass+Location.2 , data=dataframe)),#6
                      
                      AICc(lm(y~allSmolt*HakeBiomass+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HarborSeal+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Chinook*Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HakeBiomass+HarborSeal+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HakeBiomass+Chum+Location.2 , data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HakeBiomass+allSmolt+HarborSeal+Location.2 , data=dataframe))#2
                      
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
                   #"5. Chinook, Hatch", 
                   "6. Hatch, Hake", "7. Chinook, Hake, Hatch", "8. Herring, Hatch, Hake", "9. Herring, Hatch, Chinook", 
                   "11. Herring Hake Chinook",
                   
                   
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook", "Harbor Seal, Hatch", "Hake","1. Harbor Seal, Herring, Chinook",
                   "2.Harbor Seal,  Herring, Hake", "3.Harbor Seal,  Herring, Hatch", "4.Harbor Seal,  Chinook, Hake","5. Harbor Seal, Chinook, Hatch", 
                   "6. Harbor Seal, Hatch, Hake",
                   
                   
                   "Int", "Int HS herring", "Int harbor seal", "Int smolt chin", "Int Chum", "Int smolt seal")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight")
  return(aic.output)
}
model.selectionPREY <- model.selection.PREY2(dataPrey, n, dataPrey$TP)
length(model.selectionPREY[,1])
model.selectionPREY
x<-data.frame(model.selectionPREY)
subset(x, delAICc<=2)
PREY<-lm(TP~allSmolt+HakeBiomass , data=dataPrey)
summary(PREY)
vif(PREY)
subset(dataPrey, AA=="GLU")


########################     Plots#######################
library(dotwhisker)
library(broom)
library(dplyr)
library(colorspace)
library(ggpubr)

color<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=4)
color2<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=4)
color3<- rep(c('black','#CCA65A','#7EBA68','#6FB1E7','#D494E1'), each=3)
color4<- rep(c('black','#CCA65A','#7EBA68','#6FB1E7','#D494E1'), times=3)
color5<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=2)
color6<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), times=2)
color7<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), each=5)
color8<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1'), times=5)



ENV <- tidy(modelENV)
x <- data.frame("term" = c('Intercept',   'Location','MEI'), "estimate" =   c(ENV[1:3,2]), "std.error" = c(ENV[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
ENV.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[1,1],coef(modelENV)$AA[1,2]),
                "std.error" = c(ranef(modelENV)$AA[1,1],  ranef(modelENV)$AA[1,2]), "group" = 'random', model="Ala")
ENV.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[3,1],coef(modelENV)$AA[3,2]) ,
                "std.error" = c(ranef(modelENV)$AA[3,1], ranef(modelENV)$AA[3,2]), "group" = 'random', model="Pro")
ENV.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[2,1],coef(modelENV)$AA[2,2]), 
                "std.error" = c(ranef(modelENV)$AA[2,1], ranef(modelENV)$AA[2,2]), "group" = 'random', model="Glu")
ENV.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelENV)$AA[4,1],coef(modelENV)$AA[4,2]), 
                "std.error" = c(ranef(modelENV)$AA[4,1], ranef(modelENV)$AA[4,2]), "group" = 'random', model="Val")
ENV.VAL <- as_tibble(x)

ENV.mods <- rbind(ENV.fixed, ENV.GLU, ENV.ALA,   ENV.VAL, ENV.PRO)

ENV.plot <-small_multiple(ENV.mods) +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("A. Environmental") +
  
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10)) 
ENV.plot


modelPREY <- lm(TP~Location.2+allSmolt , data=dataPrey)
PREY <- tidy(modelPREY)
x <- data.frame("term" = c('Intercept',  'Location', "Smolts"), "estimate" = fixef(modelPREY), "std.error" = c(PREY[1:3,3]), "group" = c(rep('fixed', 3)), model="fixed")
PREY.fixed <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelPREY)$AA[1,1],coef(modelPREY)$AA[1,2]), "std.error" = c(ranef(modelPREY)$AA[1,1],  ranef(modelPREY)$AA[1,2]), "group" = 'random', model="Ala")
PREY.ALA <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelPREY)$AA[2,1],coef(modelPREY)$AA[2,2]) , "std.error" = c(ranef(modelPREY)$AA[2,1], ranef(modelPREY)$AA[2,2]), "group" = 'random', model="Glu")
PREY.GLU <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelPREY)$AA[3,1],coef(modelPREY)$AA[3,2]), "std.error" = c(ranef(modelPREY)$AA[3,1], ranef(modelPREY)$AA[3,2]), "group" = 'random', model="Pro")
PREY.PRO <- as_tibble(x)
x <- data.frame("term" = c('Intercept', 'Location'), "estimate" = c(coef(modelPREY)$AA[4,1],coef(modelPREY)$AA[4,2]), "std.error" = c(ranef(modelPREY)$AA[4,1], ranef(modelPREY)$AA[4,2]), "group" = 'random', model="Val")
PREY.VAL <- as_tibble(x)
PREY.mods <- rbind(PREY.fixed, PREY.GLU, PREY.ALA,  PREY.VAL, PREY.PRO)


PREY.plot <-small_multiple(PREY.mods) +
  theme_bw()+
  geom_point (colour=color3)+
  geom_errorbar (colour=color4, width=.2,
                 position=position_dodge(.9))+
  ylab("") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("B. Food Web") +
  # scale_y_discrete("", breaks=waiver(), labels=waiver())+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=14), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text( size = 10))
PREY.plot



pdf(file="Results/Presentation Figures/PreyDRAFT.pdf", width=5, height=5)
PREY.plot
dev.off()

pdf(file="Results/Presentation Figures/EnvDRAFT.pdf", width=5, height=5)
ENV.plot
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

