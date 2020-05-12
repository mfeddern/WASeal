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

data <- read.csv("Data/Compiled/HierarchicalData.csv")


########################     GLU Climate Models       ############################
dataCLIM <- subset(data, Year>=1960&Year<=2008 & AA=="Glu")
dataCLIM <-dataCLIM %>% select(MEI, 
                               PDO,
                               NPGO,
                               WA.SST.Su, 
                               TUMI.45,
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
ModelSelection.CLIM<- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lm(y~Location.2, data=dataframe)), 
                      AIC(lm(y~PDO+Location.2, data=dataframe)),
                      AIC(lm(y~NPGO+Location.2, data=dataframe)),
                      AIC(lm(y~MEI+Location.2, data=dataframe)),
                      AIC(lm(y~UpInAn.45.Spring+Location.2, data=dataframe)),
                      AIC(lm(y~PDO+NPGO+Location.2, data=dataframe)),#1
                      AIC(lm(y~PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AIC(lm(y~PDO+MEI+Location.2, data=dataframe)),#3
                      AIC(lm(y~UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AIC(lm(y~MEI+NPGO+Location.2, data=dataframe)),#5
                      AIC(lm(y~MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      AIC(lm(y~MEI+NPGO+UpInAn.45.Spring+Location.2, data=dataframe)),#7
                      AIC(lm(y~PDO+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#8
                      AIC(lm(y~PDO+MEI+NPGO+Location.2, data=dataframe)),#9
                      AIC(lm(y~PDO+NPGO+UpInAn.45.Spring+Location.2, data=dataframe)),
                      
                      AIC(lm(y~Location.2+WA.SST.Su, data=dataframe)), 
                      AIC(lm(y~WA.SST.Su+PDO+Location.2, data=dataframe)),
                      AIC(lm(y~WA.SST.Su+NPGO+Location.2, data=dataframe)),
                      AIC(lm(y~WA.SST.Su+MEI+Location.2, data=dataframe)),
                      AIC(lm(y~WA.SST.Su+UpInAn.45.Spring+Location.2, data=dataframe)),
                      AIC(lm(y~WA.SST.Su+PDO+NPGO+Location.2, data=dataframe)),#1
                      AIC(lm(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AIC(lm(y~WA.SST.Su+PDO+MEI+Location.2, data=dataframe)),#3
                      AIC(lm(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AIC(lm(y~WA.SST.Su+MEI+NPGO+Location.2, data=dataframe)),#5
                      AIC(lm(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      
                      AIC(lm(y~Location.2+TUMI.45, data=dataframe)), 
                      AIC(lm(y~TUMI.45+PDO+Location.2, data=dataframe)),
                      AIC(lm(y~TUMI.45+NPGO+Location.2, data=dataframe)),
                      AIC(lm(y~TUMI.45+MEI+Location.2, data=dataframe)),
                      AIC(lm(y~TUMI.45+UpInAn.45.Spring+Location.2, data=dataframe)),
                      AIC(lm(y~TUMI.45+PDO+NPGO+Location.2, data=dataframe)),#1
                      AIC(lm(y~TUMI.45+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AIC(lm(y~TUMI.45+PDO+MEI+Location.2, data=dataframe)),#3
                      AIC(lm(y~TUMI.45+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AIC(lm(y~TUMI.45+MEI+NPGO+Location.2, data=dataframe)),#5
                      AIC(lm(y~TUMI.45+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      
                      AIC(lm(y~Location.2+Col.Dis.high, data=dataframe)), 
                      AIC(lm(y~Col.Dis.high+PDO+Location.2, data=dataframe)),
                      AIC(lm(y~Col.Dis.high+NPGO+Location.2, data=dataframe)),
                      AIC(lm(y~Col.Dis.high+MEI+Location.2, data=dataframe)),
                      AIC(lm(y~Col.Dis.high+UpInAn.45.Spring+Location.2, data=dataframe)),
                      AIC(lm(y~Col.Dis.high+PDO+NPGO+Location.2, data=dataframe)),#1
                      AIC(lm(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AIC(lm(y~Col.Dis.high+PDO+MEI+Location.2, data=dataframe)),#3
                      AIC(lm(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AIC(lm(y~Col.Dis.high+MEI+NPGO+Location.2, data=dataframe)),#5
                      AIC(lm(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      
                      AIC(lm(y~Col.Dis.high+WA.SST.Su+Location.2, data=dataframe)),#4
                      AIC(lm(y~Col.Dis.high+TUMI.45+Location.2, data=dataframe)),#5
                      AIC(lm(y~WA.SST.Su+TUMI.45+Location.2, data=dataframe)),#6
                      AIC(lm(y~WA.SST.Su+TUMI.45+Col.Dis.high+Location.2, data=dataframe))#6

                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "PDO", "NPGO", "MEI", "Upwelling (Sp)","1. PDO, NPGO", "2. PDO, Upwelling (Sp)", "3. PDO, MEI", "4. NPGO, Upwelling (Sp)",
                   "5. NPGO, MEI", "6. MEI, Upwelling (Sp)", "7. NPGO, Upwelling (Sp), MEI", "8. PDO, MEI, Upwelling (Sp)", "9. PDO, MEI, NPGO", 
                   "11. PDO Upwelling (Sp) NPGO",
                   
                   "WA.SST.Su, Location", "WA.SST.Su, PDO", "WA.SST.Su, NPGO", "WA.SST.Su, MEI", "WA.SST.Su, Upwelling (Sp)","1. WA.SST.Su,  PDO, NPGO", "2. WA.SST.Su, PDO, Upwelling (Sp)", 
                   "3. WA.SST.Su, PDO, MEI", "4. WA.SST.Su, NPGO, Upwelling (Sp)","5. WA.SST.Su, NPGO, MEI", "6. WA.SST.Su, MEI, Upwelling (Sp)",
                   
                   "TUMI.45, Location", "TUMI.45, PDO", "TUMI.45, NPGO", "TUMI.45, MEI", "TUMI.45, Upwelling (Sp)","1. TUMI.45, PDO, NPGO", "2.TUMI.45,  PDO, Upwelling (Sp)",
                   "3. TUMI.45, PDO, MEI", "4.TUMI.45,  NPGO, Upwelling (Sp)","5.TUMI.45,  NPGO, MEI", "6. TUMI.45, MEI, Upwelling (Sp)",
                   
                   "Col.Dis.high, Location", "Col.Dis.high, PDO", "Col.Dis.high, NPGO", "Col.Dis.high, MEI", "Upwelling (Sp)","1. Col.Dis.high, PDO, NPGO",
                   "2.Col.Dis.high,  PDO, Upwelling (Sp)", "3.Col.Dis.high,  PDO, MEI", "4.Col.Dis.high,  NPGO, Upwelling (Sp)","5. Col.Dis.high, NPGO, MEI", 
                   "6. Col.Dis.high, MEI, Upwelling (Sp)",
                   
                   "Col.Dis.high, WA.SST.Su", " Col.Dis.high TUMI.45", "WA.SST.Su, TUMI.45", "WA.SST.Su, TUMI.45, Col.Dis.high")
  
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
subset(x, delAICc<=1.96)

CLIM.mod.2 <- lm(TP~PDO+Location.2, data=dataCLIM)
CLIM.mod.6 <- lm(TP~PDO+UpInAn.45.Spring+Location.2, data=dataCLIM)
CLIM.mod.4 <- lm(TP~PDO+WA.SST.Su+Location.2, data=dataCLIM)
CLIM.mod.1 <- lm(TP~PDO+TUMI.45+Location.2, data=dataCLIM)
CLIM.mod.5 <- lm(TP~PDO+TUMI.45+NPGO+Location.2, data=dataCLIM)
CLIM.mod.3 <- lm(TP~PDO+TUMI.45+MEI+Location.2, data=dataCLIM)
CLIM.mod.7 <- lm(TP~PDO+UpInAn.45.Spring+Col.Dis.high+Location.2, data=dataCLIM)

summary(lm(TP~PDO+NPGO+TUMI.45+Location.2, data=dataCLIM))
vif(lm(TP~PDO+NPGO+TUMI.45+Location.2, data=dataCLIM))
vif(CLIM.mod.1)

legend.size <- 12
CLIM <- list(CLIM.mod.1, CLIM.mod.2, CLIM.mod.3, CLIM.mod.4,
             CLIM.mod.5, CLIM.mod.6, CLIM.mod.7)
CLIM.aic <- c(AICc(CLIM.mod.1), AICc(CLIM.mod.2), AICc(CLIM.mod.3),AICc(CLIM.mod.4),
                 AICc(CLIM.mod.5), AICc(CLIM.mod.6),AICc(CLIM.mod.7))
w.CLIM <- round(akaike.weights(CLIM.aic)$weights, digits=2)
phe.CLIM.Best <-  CLIM.mod.1

CLIM.plot <-dwplot(CLIM, 
                      vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Col.Dis.high = "Columbia River Discharge (High)",                       
                       Location.2Inland = "Subregion (Salish Sea)",
                       UpInAn.45.Spring="Spring Upwelling",
                       WA.SST.Su="SST Summer",
                       TUMI.45="TUMI")) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Trophic Position", subtitle="Environmental") +
  scale_color_manual(name="Model Weights", labels = w.CLIM, values=hcl.colors(7, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 2, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))


########################    GLU PREY Models       ############################

dataPrey <- subset(data, Year>=1973&Year<=2008 & AA=="Glu")
#dataPrey <- subset(data, Year>=1973&Year<=2008 & AA=="VAL")
dataPrey <-dataPrey %>% select(allSmolt, 
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

ModelSelection.PREY <- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lm(y~Location.2, data=dataframe)), 
                      AIC(lm(y~Herring.Biomass+Location.2, data=dataframe)),
                      AIC(lm(y~Chinook+Location.2, data=dataframe)),
                      AIC(lm(y~allSmolt+Location.2, data=dataframe)),
                      AIC(lm(y~HakeBiomass+Location.2, data=dataframe)),
                      AIC(lm(y~Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AIC(lm(y~Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AIC(lm(y~Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AIC(lm(y~HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AIC(lm(y~allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AIC(lm(y~allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      AIC(lm(y~allSmolt+Chinook+HakeBiomass+Location.2, data=dataframe)),#7
                      AIC(lm(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2, data=dataframe)),#8
                      AIC(lm(y~Herring.Biomass+allSmolt+Chinook+Location.2, data=dataframe)),#9
                      AIC(lm(y~Herring.Biomass+Chinook+HakeBiomass+Location.2, data=dataframe)),
                      
                      AIC(lm(y~Location.2+Chum, data=dataframe)), 
                      AIC(lm(y~Chum+Herring.Biomass+Location.2, data=dataframe)),
                      AIC(lm(y~Chum+Chinook+Location.2, data=dataframe)),
                      AIC(lm(y~Chum+allSmolt+Location.2, data=dataframe)),
                      AIC(lm(y~Chum+HakeBiomass+Location.2, data=dataframe)),
                      AIC(lm(y~Chum+Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AIC(lm(y~Chum+Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AIC(lm(y~Chum+Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AIC(lm(y~Chum+HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AIC(lm(y~Chum+allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AIC(lm(y~Chum+allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      
                      AIC(lm(y~Location.2+Coho, data=dataframe)), 
                     AIC(lm(y~Coho+Herring.Biomass+Location.2, data=dataframe)),
                      AIC(lm(y~Coho+Chinook+Location.2, data=dataframe)),
                      AIC(lm(y~Coho+allSmolt+Location.2, data=dataframe)),
                      AIC(lm(y~Coho+HakeBiomass+Location.2, data=dataframe)),
                      AIC(lm(y~Coho+Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AIC(lm(y~Coho+Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AIC(lm(y~Coho+Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AIC(lm(y~Coho+HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AIC(lm(y~Coho+allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AIC(lm(y~Coho+allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      
                      AIC(lm(y~Location.2+HarborSeal, data=dataframe)), 
                      AIC(lm(y~HarborSeal+Herring.Biomass+Location.2, data=dataframe)),
                      AIC(lm(y~HarborSeal+Chinook+Location.2, data=dataframe)),
                      AIC(lm(y~HarborSeal+allSmolt+Location.2, data=dataframe)),
                      AIC(lm(y~HarborSeal+HakeBiomass+Location.2, data=dataframe)),
                      AIC(lm(y~HarborSeal+Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AIC(lm(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AIC(lm(y~HarborSeal+Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AIC(lm(y~HarborSeal+HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AIC(lm(y~HarborSeal+allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AIC(lm(y~HarborSeal+allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      
                      AIC(lm(y~HarborSeal+Chum+Location.2, data=dataframe)),#4
                      AIC(lm(y~HarborSeal+Coho+Location.2, data=dataframe)),#5
                      AIC(lm(y~Chum+Coho+Location.2, data=dataframe)),#6
                      AIC(lm(y~Chum+Coho+HarborSeal+Location.2, data=dataframe))#6
                      
                      
                      
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
n<- 16

model.selectionPREY <- ModelSelection.PREY(dataPrey, n, dataPrey$TP)
x<-data.frame(model.selectionPREY)
subset(x, delAICc<=1.95)

summary(lm(TP~Herring.Biomass+Location.2, data=dataPrey))
summary(lm(TP~HarborSeal+Herring.Biomass+Chinook+Location.2, data=dataPrey))
vif(lm(TP~HarborSeal+Herring.Biomass+Chinook+Location.2, data=dataPrey))
head(dataPrey)


PREY.mod.1 <- lm(TP~Herring.Biomass+Location.2, data=dataPrey)
PREY.mod.4 <- lm(TP~HakeBiomass+Location.2, data=dataPrey)
PREY.mod.5 <- lm(TP~Herring.Biomass+Chinook+Location.2, data=dataPrey)
PREY.mod.6 <- lm(TP~Herring.Biomass+allSmolt+Location.2, data=dataPrey)
PREY.mod.2 <- lm(TP~Herring.Biomass+Coho+Location.2, data=dataPrey)
PREY.mod.3 <- lm(TP~Herring.Biomass+HarborSeal+Location.2, data=dataPrey)

legend.size <- 12
PREY <- list(PREY.mod.1, PREY.mod.2, PREY.mod.3, PREY.mod.4,
             PREY.mod.5, PREY.mod.6)
PREY.aic <- c(AICc(PREY.mod.1), AICc(PREY.mod.2), AICc(PREY.mod.3),AICc(PREY.mod.4),
              AICc(PREY.mod.5), AICc(PREY.mod.6))
w.PREY <- round(akaike.weights(PREY.aic)$weights, digits=2)
phe.PREY.Best <-  PREY.mod.1

PREY.plot <-dwplot(PREY, 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Herring.Biomass = "Herring",                       
                       Location.2Inland = "Subregion (Salish Sea)",
                       HarborSeal="Harbor Seal",
                       allSmolt="Smolts",
                       TUMI.45="TUMI")) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Trophic Position", subtitle="Food Web") +
  scale_color_manual(name="Model Weights", labels = w.PREY, values=hcl.colors(7, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 2, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))



########################     GLU Nutrient Models       ############################
dataNutrient <- subset(data, Year>=1928&Year<=2008 & AA=="Glu")
dataNut <-dataNutrient %>% select(PHE.mean,
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

ModelSelection.WA3 <- function(dataframe,n, y) {
  
  aic.output <- rbind(
    AIC(lm(y~Location.2, data=dataframe)),
    AIC(lm(y~PHE.norm+d13C.norm+Location.2, data=dataframe)),
    AIC(lm(y~PHE.norm+Location.2, data=dataframe)),
    AIC(lm(y~d13C.norm+Location.2, data=dataframe))
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
model.selectionNUTRIENT <- ModelSelection.WA3(dataNut, n, dataNut$TP.norm)
model.selectionNUTRIENT
summary(lm(TP~PHE.norm+d13C.norm+Location.2, data=dataNut))
vif(lm(TP~PHE.norm+d13C.norm+Location.2, data=dataNut))


NUTR.mod.2 <- lm(TP~d13C.s+PHE.mean+Location.2, data=dataNut)
NUTR.mod.1<- lm(TP~PHE.mean+Location.2, data=dataNut)

legend.size <- 12
NUTR <- list(NUTR.mod.1, NUTR.mod.2)
NUTR.aic <- c(AICc(NUTR.mod.1), AICc(NUTR.mod.2))
w.NUTR <- round(akaike.weights(NUTR.aic)$weights, digits=2)
NUTR.Best <-  NUTR.mod.1

NUTR.plot <-dwplot(NUTR, 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(#PHE.mean = expression(paste(delta^15, "N Phenylalanine")),                       
                       Location.2Inland = "Subregion (Salish Sea)"))+
                       #d13C.s=expression(paste(delta^13, "C")))) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Trophic Position", subtitle="Nutrient") +
  scale_color_manual(name="Model Weights", labels = w.NUTR, values=hcl.colors(3, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 2, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))



###############Compile plots  ##############

pdf(file="Results/Figures/CoefPlot.Glu.pdf", width=11, height=12)
ggarrange(CLIM.plot, NUTR.plot,PREY.plot + rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 1, nrow = 3, align="hv")

dev.off()


