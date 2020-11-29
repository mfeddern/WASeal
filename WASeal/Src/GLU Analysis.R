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

data <- read.csv("Data/Compiled/HierarchicalFiles/HierarchicalData_beta1.csv")


########################     GLU Climate Models       ############################
dataCLIM <- subset(data, 
                   AA=="Glu")
dataCLIM <-dataCLIM %>% select(MEI, 
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
ModelSelection.CLIM<- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lm(y~Location.2, data=dataframe)), 
                      AICc(lm(y~PDO+Location.2, data=dataframe)),
                      AICc(lm(y~NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~MEI+Location.2, data=dataframe)),
                      AICc(lm(y~UpInAn.45.Spring+Location.2, data=dataframe)),
                      AICc(lm(y~PDO+NPGO+Location.2, data=dataframe)),#1
                      AICc(lm(y~PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AICc(lm(y~PDO+MEI+Location.2, data=dataframe)),#3
                      AICc(lm(y~UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AICc(lm(y~MEI+NPGO+Location.2, data=dataframe)),#5
                      AICc(lm(y~MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      AICc(lm(y~MEI+NPGO+UpInAn.45.Spring+Location.2, data=dataframe)),#7
                      AICc(lm(y~PDO+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#8
                      AICc(lm(y~PDO+MEI+NPGO+Location.2, data=dataframe)),#9
                      AICc(lm(y~PDO+NPGO+UpInAn.45.Spring+Location.2, data=dataframe)),
                      
                      AICc(lm(y~Location.2+WA.SST.Su, data=dataframe)), 
                      AICc(lm(y~WA.SST.Su+PDO+Location.2, data=dataframe)),
                      AICc(lm(y~WA.SST.Su+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~WA.SST.Su+MEI+Location.2, data=dataframe)),
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Spring+Location.2, data=dataframe)),
                      AICc(lm(y~WA.SST.Su+PDO+NPGO+Location.2, data=dataframe)),#1
                      AICc(lm(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AICc(lm(y~WA.SST.Su+PDO+MEI+Location.2, data=dataframe)),#3
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AICc(lm(y~WA.SST.Su+MEI+NPGO+Location.2, data=dataframe)),#5
                      AICc(lm(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      
                      AICc(lm(y~Location.2+UpInAn.45.Summer, data=dataframe)), 
                      AICc(lm(y~UpInAn.45.Summer+PDO+Location.2, data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+MEI+Location.2, data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+UpInAn.45.Spring+Location.2, data=dataframe)),
                      AICc(lm(y~UpInAn.45.Summer+PDO+NPGO+Location.2, data=dataframe)),#1
                      AICc(lm(y~UpInAn.45.Summer+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AICc(lm(y~UpInAn.45.Summer+PDO+MEI+Location.2, data=dataframe)),#3
                      AICc(lm(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AICc(lm(y~UpInAn.45.Summer+MEI+NPGO+Location.2, data=dataframe)),#5
                      AICc(lm(y~UpInAn.45.Summer+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      
                      AICc(lm(y~Location.2+Col.Dis.high, data=dataframe)), 
                      AICc(lm(y~Col.Dis.high+PDO+Location.2, data=dataframe)),
                      AICc(lm(y~Col.Dis.high+NPGO+Location.2, data=dataframe)),
                      AICc(lm(y~Col.Dis.high+MEI+Location.2, data=dataframe)),
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Spring+Location.2, data=dataframe)),
                      AICc(lm(y~Col.Dis.high+PDO+NPGO+Location.2, data=dataframe)),#1
                      AICc(lm(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2, data=dataframe)),#2
                      AICc(lm(y~Col.Dis.high+PDO+MEI+Location.2, data=dataframe)),#3
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2, data=dataframe)),#4
                      AICc(lm(y~Col.Dis.high+MEI+NPGO+Location.2, data=dataframe)),#5
                      AICc(lm(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2, data=dataframe)),#6
                      
                      AICc(lm(y~Col.Dis.high+WA.SST.Su+Location.2, data=dataframe)),#4
                      AICc(lm(y~Col.Dis.high+UpInAn.45.Summer+Location.2, data=dataframe)),#5
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Summer+Location.2, data=dataframe)),#6
                      AICc(lm(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2, data=dataframe))#6

                      
                      
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

ModelSelection.CLIM.hier<- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|Location.2), data=dataframe)), 
                      AICc(lmer(y~PDO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~NPGO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~MEI+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Spring+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~PDO+NPGO+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~PDO+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~PDO+MEI+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~UpInAn.45.Spring+NPGO+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~MEI+NPGO+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~MEI+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#6
                      AICc(lmer(y~MEI+NPGO+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#7
                      AICc(lmer(y~PDO+MEI+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#8
                      AICc(lmer(y~PDO+MEI+NPGO+(1|Location.2), data=dataframe)),#9
                      AICc(lmer(y~PDO+NPGO+UpInAn.45.Spring+(1|Location.2), data=dataframe)),
                      
                      AICc(lmer(y~(1|Location.2)+WA.SST.Su, data=dataframe)), 
                      AICc(lmer(y~WA.SST.Su+PDO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~WA.SST.Su+NPGO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~WA.SST.Su+MEI+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~WA.SST.Su+PDO+NPGO+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~WA.SST.Su+PDO+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~WA.SST.Su+PDO+MEI+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+NPGO+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~WA.SST.Su+MEI+NPGO+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~WA.SST.Su+MEI+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#6
                      
                      AICc(lmer(y~(1|Location.2)+UpInAn.45.Summer, data=dataframe)), 
                      AICc(lmer(y~UpInAn.45.Summer+PDO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Summer+NPGO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Summer+MEI+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~UpInAn.45.Summer+PDO+NPGO+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~UpInAn.45.Summer+PDO+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~UpInAn.45.Summer+PDO+MEI+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~UpInAn.45.Summer+MEI+NPGO+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~UpInAn.45.Summer+MEI+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#6
                      
                      AICc(lmer(y~(1|Location.2)+Col.Dis.high, data=dataframe)), 
                      AICc(lmer(y~Col.Dis.high+PDO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+NPGO+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+MEI+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Col.Dis.high+PDO+NPGO+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~Col.Dis.high+PDO+MEI+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~Col.Dis.high+MEI+NPGO+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+(1|Location.2), data=dataframe)),#6
                      
                      AICc(lmer(y~Col.Dis.high+WA.SST.Su+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Summer+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+(1|Location.2), data=dataframe)),#6
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+(1|Location.2), data=dataframe))#6
                      
                      
                      
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
subset(x, delAICc<=1.96)

summary(lm(TP~PDO+Location.2, data=dataCLIM))


CLIM.mod.1 <- lm(TP~PDO+Location.2, data=dataCLIM)
CLIM.mod.2 <- lm(TP~PDO+Col.Dis.high+Location.2, data=dataCLIM)


legend.size <- 12
CLIM <- list(CLIM.mod.1, CLIM.mod.2)
CLIM.aic <- c(AICc(CLIM.mod.1), AICc(CLIM.mod.2))
w.CLIM <- round(akaike.weights(CLIM.aic)$weights, digits=2)
phe.CLIM.Best <-  CLIM.mod.1

summary(CLIM.mod.1)
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
  scale_color_manual(name="Model Weights", labels = w.CLIM, values=hcl.colors(3, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0.5, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))


########################    GLU PREY Models       ############################

dataPrey <- subset(data, AA=="Glu")
#dataPrey <- subset(data, Year>=1973&Year<=2008 )
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

n<- 16
ModelSelection.PREY <- function(dataframe,n,y) {
  
  aic.output <- rbind(AICc(lm(y~Location.2, data=dataframe)), 
                      AICc(lm(y~Herring.Biomass+Location.2, data=dataframe)),
                      AICc(lm(y~Chinook+Location.2, data=dataframe)),
                      AICc(lm(y~allSmolt+Location.2, data=dataframe)),
                      AICc(lm(y~HakeBiomass+Location.2, data=dataframe)),
                      AICc(lm(y~Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AICc(lm(y~Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AICc(lm(y~Herring.Biomass*HakeBiomass+Location.2, data=dataframe)),#2
                      
                      AICc(lm(y~Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AICc(lm(y~HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AICc(lm(y~allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AICc(lm(y~allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      AICc(lm(y~allSmolt+Chinook+HakeBiomass+Location.2, data=dataframe)),#7
                      AICc(lm(y~Herring.Biomass+HakeBiomass+allSmolt+Location.2, data=dataframe)),#8
                      AICc(lm(y~Herring.Biomass+allSmolt+Chinook+Location.2, data=dataframe)),#9
                      AICc(lm(y~Herring.Biomass+HakeBiomass+Chinook+Location.2, data=dataframe)),
                      
                      AICc(lm(y~Location.2+Chum, data=dataframe)), 
                      AICc(lm(y~Chum+Herring.Biomass+Location.2, data=dataframe)),
                      AICc(lm(y~Chum+Chinook+Location.2, data=dataframe)),
                      AICc(lm(y~Chum+allSmolt+Location.2, data=dataframe)),
                      AICc(lm(y~Chum+HakeBiomass+Location.2, data=dataframe)),
                      AICc(lm(y~Chum+Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AICc(lm(y~Chum+Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AICc(lm(y~Chum+Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AICc(lm(y~Chum+HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AICc(lm(y~Chum+allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AICc(lm(y~Chum+allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      
                      AICc(lm(y~Location.2+Coho, data=dataframe)), 
                      AICc(lm(y~Coho+Herring.Biomass+Location.2, data=dataframe)),
                      AICc(lm(y~Coho+Chinook+Location.2, data=dataframe)),
                      AICc(lm(y~Coho+allSmolt+Location.2, data=dataframe)),
                      AICc(lm(y~Coho+HakeBiomass+Location.2, data=dataframe)),
                      AICc(lm(y~Coho+Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AICc(lm(y~Coho+Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AICc(lm(y~Coho+Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AICc(lm(y~Coho+HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AICc(lm(y~Coho+allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AICc(lm(y~Coho+allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      
                      AICc(lm(y~Location.2+HarborSeal, data=dataframe)), 
                      AICc(lm(y~HarborSeal+Herring.Biomass+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+Chinook+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+allSmolt+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+HakeBiomass+Location.2, data=dataframe)),
                      AICc(lm(y~HarborSeal+Herring.Biomass+Chinook+Location.2, data=dataframe)),#1
                      AICc(lm(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2, data=dataframe)),#2
                      AICc(lm(y~HarborSeal+Herring.Biomass+allSmolt+Location.2, data=dataframe)),#3
                      AICc(lm(y~HarborSeal+HakeBiomass+Chinook+Location.2, data=dataframe)),#4
                      AICc(lm(y~HarborSeal+allSmolt+Chinook+Location.2, data=dataframe)),#5
                      AICc(lm(y~HarborSeal+allSmolt+HakeBiomass+Location.2, data=dataframe)),#6
                      
                      AICc(lm(y~HarborSeal+Chum+Location.2, data=dataframe)),#4
                      AICc(lm(y~HarborSeal+Coho+Location.2, data=dataframe)),#5
                      AICc(lm(y~Chum+Coho+Location.2, data=dataframe)),#6
                      AICc(lm(y~Chum+Coho+HarborSeal+Location.2, data=dataframe))#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location",
                   "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "Herring*Hake","3. Herring, Hatch", "4. Chinook, Hake",
                   "5. Chinook, Hatch", "6. Hatch, Hake", "7. Chinook, Hake, Hatch", "8. Herring, Hatch, Hake", "9. Herring, Hatch, Chinook", 
                   "11. Herring Hake Chinook",
                   
                   "Chum, Location", "Chum, Herring", "Chum, Chinook", "Chum, Hatch", "Chum, Hake","1. Chum,  Herring, Chinook", "2. Chum, Herring, Hake", 
                   "3. Chum, Herring, Hatch", "4. Chum, Chinook, Hake","5. Chum, Chinook, Hatch", "6. Chum, Hatch, Hake",
                   
                   "Coho, Location", "Coho, Herring", "Coho, Chinook", "Coho, Hatch", "Coho, Hake","1. Coho, Herring, Chinook", "2.Coho,  Herring, Hake",
                   "3. Coho, Herring, Hatch", "4.Coho,  Chinook, Hake","5.Coho,  Chinook, Hatch", "6. Coho, Hatch, Hake",
                   
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook", "Harbor Seal, Hatch", "Harbor Seal Hake","Harbor Seal, Herring, Chinook",
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

ModelSelection.PREYhier <- function(dataframe,n,y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|Location.2), data=dataframe)), 
                      AICc(lmer(y~Herring.Biomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Chinook+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~allSmolt+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~HakeBiomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Herring.Biomass+Chinook+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~Herring.Biomass*HakeBiomass+(1|Location.2), data=dataframe)),#2
                      
                      AICc(lmer(y~Herring.Biomass+allSmolt+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~HakeBiomass+Chinook+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~allSmolt+Chinook+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~allSmolt+HakeBiomass+(1|Location.2), data=dataframe)),#6
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+(1|Location.2), data=dataframe)),#7
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+allSmolt+(1|Location.2), data=dataframe)),#8
                      AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+(1|Location.2), data=dataframe)),#9
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+Chinook+(1|Location.2), data=dataframe)),
                      
                      AICc(lmer(y~(1|Location.2)+Chum, data=dataframe)), 
                      AICc(lmer(y~Chum+Herring.Biomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Chum+Chinook+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Chum+allSmolt+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Chum+HakeBiomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Chum+Herring.Biomass+Chinook+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~Chum+Herring.Biomass+HakeBiomass+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~Chum+Herring.Biomass+allSmolt+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~Chum+HakeBiomass+Chinook+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~Chum+allSmolt+Chinook+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~Chum+allSmolt+HakeBiomass+(1|Location.2), data=dataframe)),#6
                      
                      AICc(lmer(y~(1|Location.2)+Coho, data=dataframe)), 
                      AICc(lmer(y~Coho+Herring.Biomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Coho+Chinook+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Coho+allSmolt+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Coho+HakeBiomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~Coho+Herring.Biomass+Chinook+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~Coho+Herring.Biomass+HakeBiomass+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~Coho+Herring.Biomass+allSmolt+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~Coho+HakeBiomass+Chinook+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~Coho+allSmolt+Chinook+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~Coho+allSmolt+HakeBiomass+(1|Location.2), data=dataframe)),#6
                      
                      AICc(lmer(y~(1|Location.2)+HarborSeal, data=dataframe)), 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Chinook+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~HarborSeal+allSmolt+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~HarborSeal+HakeBiomass+(1|Location.2), data=dataframe)),
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+(1|Location.2), data=dataframe)),#1
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+(1|Location.2), data=dataframe)),#2
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+(1|Location.2), data=dataframe)),#3
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+allSmolt+Chinook+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+(1|Location.2), data=dataframe)),#6
                      
                      AICc(lmer(y~HarborSeal+Chum+(1|Location.2), data=dataframe)),#4
                      AICc(lmer(y~HarborSeal+Coho+(1|Location.2), data=dataframe)),#5
                      AICc(lmer(y~Chum+Coho+(1|Location.2), data=dataframe)),#6
                      AICc(lmer(y~Chum+Coho+HarborSeal+(1|Location.2), data=dataframe))#6
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location",
                   "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "Herring*Hake","3. Herring, Hatch", "4. Chinook, Hake",
                   "5. Chinook, Hatch", "6. Hatch, Hake", "7. Chinook, Hake, Hatch", "8. Herring, Hatch, Hake", "9. Herring, Hatch, Chinook", 
                   "11. Herring Hake Chinook",
                   
                   "Chum, Location", "Chum, Herring", "Chum, Chinook", "Chum, Hatch", "Chum, Hake","1. Chum,  Herring, Chinook", "2. Chum, Herring, Hake", 
                   "3. Chum, Herring, Hatch", "4. Chum, Chinook, Hake","5. Chum, Chinook, Hatch", "6. Chum, Hatch, Hake",
                   
                   "Coho, Location", "Coho, Herring", "Coho, Chinook", "Coho, Hatch", "Coho, Hake","1. Coho, Herring, Chinook", "2.Coho,  Herring, Hake",
                   "3. Coho, Herring, Hatch", "4.Coho,  Chinook, Hake","5.Coho,  Chinook, Hatch", "6. Coho, Hatch, Hake",
                   
                   "Harbor Seal, Location", "Harbor Seal, Herring", "Harbor Seal, Chinook", "Harbor Seal, Hatch", "Harbor Seal Hake","Harbor Seal, Herring, Chinook",
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

model.selectionPREY <- ModelSelection.PREYhier(dataPrey, n, dataPrey$TP)
x<-data.frame(model.selectionPREY)
subset(x, delAICc<=1.97)




summary(lm(TP~HakeBiomass+Herring.Biomass+Location.2, data=dataPrey))
vif(lm(TP~HarborSeal+Herring.Biomass+Chinook+Location.2, data=dataPrey))
head(dataPrey)

dataPrey$Year

PREY.mod.1 <-lm(TP~HakeBiomass+Location.2, data=dataPrey)
PREY.mod.2 <-lm(TP~HakeBiomass+Herring.Biomass+Location.2, data=dataPrey)
summary(PREY.mod.2)


PREY <- list(PREY.mod.1)
PREY.aic <- c(AICc(PREY.mod.1))
              
w.PREY <- round(akaike.weights(PREY.aic)$weights, digits=2)
phe.PREY.Best <-  PREY.mod.1
col<-c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1')
legend.size<-12


 INT.Plot<- sjPlot::plot_model(PREY.mod.1, type = "int", col=c(col[1], col[2]))+
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



PREY.plot <-dwplot(PREY, 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Herring.Biomass = "Herring",                       
                       Location.2Inland = "Subregion (Salish Sea)",
                       HarborSeal="Harbor Seal",
                       allSmolt="Smolts",
                       HakeBiomass="Hake")) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Trophic Position", subtitle="Food Web") +
  scale_color_manual(name="Model Weights", labels = w.PREY, values=hcl.colors(7, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0.5, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))



########################     GLU Nutrient Models       ############################
dataNutrient <- subset(data,  AA=="Glu")
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
    AICc(lm(y~Location.2, data=dataframe)),
    AICc(lm(y~PHE.norm+d13C.norm+Location.2, data=dataframe)),
    AICc(lm(y~PHE.norm+Location.2, data=dataframe)),
    AICc(lm(y~d13C.norm+Location.2, data=dataframe))
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

vif(lm(TP~PHE.norm+d13C.norm+Location.2, data=dataNut))


#NUTR.mod.2 <- lm(TP~d13C.s+PHE.mean+Location.2, data=dataNut)
NUTR.mod.1<- lm(TP~PHE.mean+Location.2, data=dataNut)

legend.size <- 12
NUTR <- list(NUTR.mod.1)
NUTR.aic <- c(AICc(NUTR.mod.1))
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
        plot.margin = margin(0.25, 0.5, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))


###############Save plots  ##############

pdf(file="Results/Presentation Figures/GLUNutrientDRAFT.pdf", width=5, height=5)
NUTR.plot
dev.off()

pdf(file="Results/Presentation Figures/PreyDRAFT.pdf", width=5, height=5)
PREY.plot
dev.off()

pdf(file="Results/Presentation Figures/EnvDRAFT.pdf", width=5, height=5)
CLIM.plot
dev.off()

pdf(file="Results/Presentation Figures/FullDRAFT.pdf", width=5, height=5)
INT.Plot
dev.off()


###############Compile plots  ##############

pdf(file="Results/Figures/CoefPlot.Glu2.pdf", width=14, height=8)
ggarrange(CLIM.plot, PREY.plot,NUTR.plot, INT.Plot + rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 2, align="hv",
          heights =c(5,5))


dev.off()


