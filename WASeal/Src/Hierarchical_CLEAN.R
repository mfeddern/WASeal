
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
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe, REML=F)), #1
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe, REML=F)), #2
                      AICc(lmer(y~PDO+Location.2+(1|AA), data=dataframe, REML=F)),#3
                      AICc(lmer(y~NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#4
                      AICc(lmer(y~MEI+Location.2+(1|AA), data=dataframe, REML=F)),#5
                      AICc(lmer(y~UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#6
                      AICc(lmer(y~PDO+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#7
                      AICc(lmer(y~PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#8
                      AICc(lmer(y~UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#9
                      AICc(lmer(y~MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#10
                      AICc(lmer(y~Location.2+(1|AA)+WA.SST.Su, data=dataframe, REML=F)), #11
                      AICc(lmer(y~WA.SST.Su+PDO+Location.2+(1|AA), data=dataframe, REML=F)),#12
                      AICc(lmer(y~WA.SST.Su+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#13
                      AICc(lmer(y~WA.SST.Su+MEI+Location.2+(1|AA), data=dataframe, REML=F)),#14
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#15
                      AICc(lmer(y~WA.SST.Su+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#16
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#17
                      AICc(lmer(y~WA.SST.Su+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#18
                      AICc(lmer(y~Location.2+(1|AA)+UpInAn.45.Summer, data=dataframe, REML=F)), #19
                      AICc(lmer(y~UpInAn.45.Summer+PDO+Location.2+(1|AA), data=dataframe, REML=F)),#20
                      AICc(lmer(y~UpInAn.45.Summer+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#21
                      AICc(lmer(y~UpInAn.45.Summer+MEI+Location.2+(1|AA), data=dataframe, REML=F)),#22
                      AICc(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#23
                      AICc(lmer(y~Location.2+(1|AA)+Col.Dis.high, data=dataframe, REML=F)), #24
                      AICc(lmer(y~Col.Dis.high+PDO+Location.2+(1|AA), data=dataframe, REML=F)),#25
                      AICc(lmer(y~Col.Dis.high+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#26
                      AICc(lmer(y~Col.Dis.high+MEI+Location.2+(1|AA), data=dataframe, REML=F)),#27
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#28
                      AICc(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#29
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+Location.2+(1|AA), data=dataframe, REML=F)),#30
                      AICc(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+Location.2+(1|AA), data=dataframe, REML=F)),#31
                      AICc(lmer(y~Col.Dis.high+WA.SST.Su+Location.2+(1|AA), data=dataframe, REML=F)),#32
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe, REML=F)),#33
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Location.2+(1|AA), data=dataframe, REML=F)),#34
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+Location.2+(1|AA), data=dataframe, REML=F))#35

                      
                      
                      
  )
  
  model.names <- c("1. Null",
                   "2. Location Only", 
                   "3. PDO, Location", 
                   "4. NPGO, Location", 
                   "5. MEI, Location", 
                   "6. Upwelling (Spring), Location",
                   "7. NPGO, PDO, Location", 
                   "8. PDO, Upwelling (Spring), Location",  
                   "9. NPGO, Upwelling (Spring), Location",
                   "10. MEI, Upwelling (Spring), Location",  
                   "11. SST (Summer), Location", 
                   "12. SST (Summer), PDO, Location", 
                   "13. SST (Summer), NPGO, Location",
                   "14. SST (Summer), MEI, Location",
                   "15. SST (Summer), Upwelling (Spring), Location",
                   "16. SST (Summer), PDO, Upwelling (Spring), Location", 
                   "17. SST (Summer), NPGO, Upwelling (Spring), Location",
                   "18. SST (Summer), MEI, Upwelling (Spring), Location",
                   
                   "19. Upwelling (Summer), Location", 
                   "20. Upwelling (Summer), PDO, Location", 
                   "21. Upwelling (Summer), NPGO, Location", 
                   "22. Upwelling (Summer), MEI, Location",
                   "23. Upwelling (Summer),  NPGO, Upwelling (Spring), Location", 
                   "24. Columbia Discharge (High), Location", 
                   "25. Columbia Discharge (High), PDO, Location", 
                   "26. Columbia Discharge (High), NPGO, Location",
                   "27. Columbia Discharge (High), MEI, Location", 
                   "28. Upwelling (Spring), Location",
                   "29. Columbia Discharge (High),  PDO, Upwelling (Spring), Location",
                   "30. Columbia Discharge High,  NPGO, Upwelling (Spring), Location",
                   "31. Columbia Discharge High, MEI, Upwelling (Spring), Location",
                   
                   "32. Columbia Discharge (High), SST (Summer), Location", 
                   "33. Columbia Discharge (High), Upwelling (Summer), Location", 
                   "34. SST (Summer), Upwelling (Summer), Location", 
                   "35. SST (Summer), Upwelling (Summer), Columbia Discharge (High), Location")
  
 # row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- round(aic.weight1/sum(aic.weight1), digits=2)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(model.names, round(aic.output, digits=2), round(delaic, digits=2), aic.weight)
  
  colnames(aic.output)<- c("Covariates","AICc", "delAICc", "AICc Weight")
  return(aic.output)
}

Enviro.Candidate <- data.frame(Covariates = clim.1$Covariates)
sjPlot::tab_df(Enviro.Candidate,
               title = "Ocean Condition Candidate Models", 
               file = "Results/Tables/CandidateOceanCondition.doc")

model.selectionCLIM.0 <- ModelSelection.CLIM(dataCLIM.0, n, dataCLIM.0$TP)
clim.0 <-data.frame(model.selectionCLIM.0)
clim.0.ordered <- clim.0[order(clim.0$AICc),]
subset(clim.0, delAICc<=5)
clim.0.ordered[1:5,]

#### CLIMATE 1 ####
model.selectionCLIM.1 <- ModelSelection.CLIM(dataCLIM.1, n, dataCLIM.1$TP)
clim.1 <-data.frame(model.selectionCLIM.1)
clim.1.ordered <- clim.1[order(clim.1$AICc),]
subset(clim.1, delAICc<=5)
clim.1.ordered[1:10,]

modelCLIM1.1<-lmer(TP~Location.2+UpInAn.45.Summer+UpInAn.45.Spring+NPGO+(1|AA), data=dataCLIM.1, REML=F)
modelCLIM1.2<-lmer(TP~Location.2+UpInAn.45.Summer+Col.Dis.high+WA.SST.Su+(1|AA), data=dataCLIM.1, REML=F)
modelCLIM1.3<-lmer(TP~Location.2+UpInAn.45.Summer+Col.Dis.high+(1|AA), data=dataCLIM.1, REML=F)
modelCLIM1.4<-lmer(TP~Location.2+UpInAn.45.Summer+(1|AA), data=dataCLIM.1, REML=F)
modelCLIM1.5<-lmer(TP~Location.2+UpInAn.45.Summer+WA.SST.Su+(1|AA), data=dataCLIM.1, REML=F)

CLIM1 <- list(modelCLIM1.1, modelCLIM1.2, modelCLIM1.3,
              modelCLIM1.4, modelCLIM1.5)
CLIM1.aic <- round(c(AICc(modelCLIM1.1), AICc(modelCLIM1.2), AICc(modelCLIM1.3),
                     AICc(modelCLIM1.4), AICc(modelCLIM1.5)), digits=2)
w.CLIM1<- clim.1.ordered$AICc.Weight[1:5]
CLIM1.Best <-  modelCLIM1.1

legend.size <- 12
CLIM1.plot <-dwplot(CLIM1, 
                    vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       WA.SST.Su="Summer \n SST",
                       UpInAn.45.Spring="Spring \n Upwelling",
                       MEI="ENSO",
                       PDO="PDO",
                       Col.Dis.high = "Columbia River \n Discharge \n (High)" ,
                       UpInAn.45.Summer="Summer \n Upwelling",
                       NPGO="NPGO"
                       )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle(aes(label="C. Physiological Delay", fontface=plain)) +
  scale_color_manual(name="Model AIC", labels = CLIM1.aic, values=hcl.colors(7, palette = "Greens", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

CLIM1.plot


sjPlot::tab_df(clim.1.ordered[1:10,],
               title = "Physiological Delay Top 10 Models", 
               file = "Results/Tables/Clim1Top5.doc")



##### CLIMATE 2#####

model.selectionCLIM.2 <- ModelSelection.CLIM(dataCLIM.2, n, dataCLIM.2$TP)
clim.2 <-data.frame(model.selectionCLIM.2)
clim.2.ordered <- clim.2[order(clim.2$AICc),]
subset(clim.2, delAICc<=1.97)
clim.2.ordered[1:10,]

sjPlot::tab_df(clim.2.ordered[1:10,],
               title = "1-Year Ecological Delay Top 10 Models", 
               file = "Results/Tables/Clim2Top5.doc")



modelCLIM2.1<-lmer(TP~Location.2+WA.SST.Su+UpInAn.45.Spring+(1|AA), data=dataCLIM.2, REML=F)
modelCLIM2.2<-lmer(TP~Location.2+WA.SST.Su+UpInAn.45.Spring+NPGO+(1|AA), data=dataCLIM.2, REML=F)
modelCLIM2.3<-lmer(TP~Location.2+WA.SST.Su+UpInAn.45.Spring+PDO+(1|AA), data=dataCLIM.2, REML=F)
modelCLIM2.4<-lmer(TP~Location.2+WA.SST.Su+UpInAn.45.Spring+MEI+(1|AA), data=dataCLIM.2, REML=F)


CLIM2 <- list(modelCLIM2.1, modelCLIM2.2, modelCLIM2.3,
              modelCLIM2.4)
CLIM2.aic <- round(c(AICc(modelCLIM2.1), AICc(modelCLIM2.2), AICc(modelCLIM2.3),
                     AICc(modelCLIM2.4)), digits=2)
w.CLIM2<- clim.2.ordered$AICc.Weight[1:4]
CLIM1.Best <-  modelCLIM1.1

legend.size <- 12
CLIM2.plot <-dwplot(CLIM2, 
                    vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       WA.SST.Su="Summer \n SST",
                       UpInAn.45.Spring="Spring \n Upwelling",
                       MEI="ENSO",
                       PDO="PDO",
                       Col.Dis.high = "Columbia River \n Discharge \n (High)" ,
                       UpInAn.45.Summer="Summer \n Upwelling",
                       NPGO="NPGO")) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("B. 1-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = CLIM2.aic, values=hcl.colors(7, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

CLIM2.plot



#### CLIMATE 3 #####



model.selectionCLIM.3 <- ModelSelection.CLIM(dataCLIM.3, n, dataCLIM.3$TP)
clim.3 <-data.frame(model.selectionCLIM.3)
subset(clim.3, delAICc<=2)
clim.3.ordered <- clim.3[order(clim.3$AICc),]
clim.3.ordered[1:10,]

sjPlot::tab_df(clim.3.ordered[1:10,],
               title = "1-Year Ecological Delay Top 10 Models", 
               file = "Results/Tables/Clim3Top5.doc")




modelCLIM3.1<-lmer(TP~Location.2+Col.Dis.high+(1|AA), data=dataCLIM.3, REML=F)
modelCLIM3.2<-lmer(TP~Location.2+Col.Dis.high+UpInAn.45.Summer+(1|AA), data=dataCLIM.3, REML=F)
modelCLIM3.3<-lmer(TP~Location.2+Col.Dis.high+NPGO+(1|AA), data=dataCLIM.3, REML=F)
modelCLIM3.4<-lmer(TP~Location.2+Col.Dis.high+PDO+(1|AA), data=dataCLIM.3, REML=F)
modelCLIM3.5<-lmer(TP~Location.2+Col.Dis.high+MEI+(1|AA), data=dataCLIM.3, REML=F)


CLIM3 <- list(modelCLIM3.1, modelCLIM3.2, modelCLIM3.3,
              modelCLIM3.4,modelCLIM3.5)
CLIM3.aic <- round(c(AICc(modelCLIM3.1), AICc(modelCLIM3.2), AICc(modelCLIM3.3),
                     AICc(modelCLIM3.4), AICc(modelCLIM3.5)), digits=2)
w.CLIM3<- clim.3.ordered$AICc.Weight[1:5]
CLIM3.Best <-  modelCLIM3.1

legend.size <- 12
CLIM3.plot <-dwplot(CLIM3, 
                    vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       WA.SST.Su="Summer \n SST",
                       UpInAn.45.Spring="Spring \n Upwelling",
                       MEI="ENSO",
                       PDO="PDO",
                       Col.Dis.high = "Columbia River \n Discharge \n (High)" ,
                       UpInAn.45.Summer="Summer \n Upwelling",
                       NPGO="NPGO")) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("A. 2-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = CLIM3.aic, values=hcl.colors(7, palette = "Magenta", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

CLIM3.plot




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

model.selection.PREY<- function(dataframe,n, y) {
  
  aic.output <- rbind(AICc(lmer(y~(1|AA), data=dataframe, REML=F)), #1, null
                      AICc(lmer(y~Location.2+(1|AA), data=dataframe, REML=F)), #2 Location Only
                      AICc(lmer(y~Herring.Biomass+Location.2+(1|AA), data=dataframe, REML=F)), #3. Herring Biomass
                      AICc(lmer(y~Chinook+Location.2+(1|AA), data=dataframe, REML=F)),#4. Chinook Escapements
                      AICc(lmer(y~allSmolt+Location.2+(1|AA), data=dataframe, REML=F)), #5. All smolts
                      AICc(lmer(y~HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)), #6. Hake Biomass
                      AICc(lmer(y~Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe, REML=F)),#7. Herring Biomass, Chinook Escapements
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+(1|AA), data=dataframe, REML=F)),#8. Herring Biomass, Hake Biomass
                      AICc(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe, REML=F)),#9. Herring Biomass, all smolt
                      AICc(lmer(y~HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe, REML=F)),#10. Hake Bimass, Chinook Escapement
                      AICc(lmer(y~allSmolt+Chinook+Location.2+(1|AA), data=dataframe, REML=F)),#11. all smolt, Chinook Escapements
                      AICc(lmer(y~allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)),#12. all smolt, Hake Biomass
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)),#13. all smolt, Chinook escapent, Hake Biomas
                      AICc(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)),#14. Herring Biomass, all smolts, Hake Biomass
                      AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|AA), data=dataframe, REML=F)),#15.Herring biomass, all smolts, Chinook escapements
                      AICc(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)), #16. Herring Biomass, Chinook escapements, Hake biomass
                      AICc(lmer(y~HarborSeal+Location.2+(1|AA), data=dataframe, REML=F)), #17. Harbor seal 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA), data=dataframe, REML=F)),#18. Harbor seal, herring biomass
                      AICc(lmer(y~HarborSeal+Chinook+Location.2+(1|AA), data=dataframe, REML=F)), #19. Harbor seal, Chinook escapement
                      AICc(lmer(y~HarborSeal+allSmolt+Location.2+(1|AA), data=dataframe, REML=F)),#20. Harbor seal, all smolt
                      AICc(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)),#21. Harbor seal, hake biomass
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe, REML=F)),#22. Harbor seal, herring biomass, Chinook escapement
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)),#23. Harbor seal, herring biomass, hakr biomass
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|AA), data=dataframe, REML=F)),#24. Harbor seal, herring biomass, all smolts
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe, REML=F)),#25. Harbor seal, hake biomass, Chinool escapement
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F))#26. Harbor seal, all smolts, hake biomass
                    #  AICc(lmer(y~Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)) #27. Hake:Herring biomass
                      
                      
                      
  )
  
  model.names <- c("1. Null",
                   "2. Location Only",
                   "3. Herring Spawning Biomass, Location", 
                   "4. Chinook Escapements, Location", 
                   "5. Chinook Smolts Production, Location", 
                   "6. Hake Spawning Biomass, Location",
                   "7. Herring Spawning Biomass, Chinook Escapements, Location", 
                   "8. Herring Spawning Biomass, Hake Spawning Biomass, Location", 
                   "9. Herring Spawning Biomass, Chinook Smolt Production, Location", 
                   "10. Chinook Escapements, Hake Spawning Biomass, Location",
                   "11. Chinook Escapments, Chinook Smolt Production, Location",
                   "12. Chinook Smolt Production, Hake Spawning Biomass, Location", 
                   "13. Chinook Escapement, Chinook Smolt Production, Hake Spawning Biomass, Location", 
                   "14. Herring Spawning Biomass, Chinook Smolt Production, Hake Spawning Biomass, Location",
                   "15. Chinook Escapements, Chinook Smolt Production, Herring Spawning Biomass, Location",
                   "16. Herring Spawning Biomass, Hake Spawning Biomass, Chinook Escapements, Location",
                   "17. Harbor Seal Abundance, Location",
                   "18. Harbor Seal Abundance, Herring Spawning Biomass, Location", 
                   "19. Harbor Seal Abundance, Chinook Escapements, Location", 
                   "20. Harbor Seal Abundance, Chinook Smolt Production, Location", 
                   "21. Harbor seal Abundance, Hake Spawning iomass, Location",
                   "22. Harbor Seal Abundance, Herring biomass, Chinook Escapements, Location",
                   "23. Harbor Seal Abundance, Herring Spawning Biomass, Hake Spawning Biomass, Location", 
                   "24. Harbor Seal Abundance, Herring Spawning Biomass, Chinook Smolt Production, Location", 
                   "25. Harbor Seal Abundance, Chinook Escapements, Hake Spawning Biomass, Location",
                   "26. Harbor Seal Abundance, Chinook Smolt Production, Hake Spawning Biomass, Location"
                  # "27. Hake Spawning Biomass, Herring Spawning Biomass, Hake and Herring Interaction, Location"
                   )
  
  #row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(model.names, round(aic.output, digits=2), round(delaic, digits=2), round(aic.weight, digits=2))
  
  colnames(aic.output)<- c("Covariates","AICc", "delAICc", "AICc Weight")
  return(aic.output)
}





#######PREY LAG 1 #######
Prey.Candidate <- data.frame(Covariates = prey1$Covariates)
sjPlot::tab_df(Prey.Candidate,
               title = "Food Web Candidate Models", 
               file = "Results/Tables/CandidatePrey.doc")


model.selectionPREY.1 <- model.selection.PREY(dataPrey.1, n, dataPrey.1$TP)
prey1<-data.frame(model.selectionPREY.1)
prey1.ordered <- prey1[order(prey1$AICc),]
prey1.ordered[1:10,]
modelPREY1.1<-lmer(TP~Location.2+allSmolt+Herring.Biomass+HakeBiomass+(1|AA), data=dataPrey.1, REML=F)
modelPREY1.2<-lmer(TP~Location.2+allSmolt+HakeBiomass+(1|AA), data=dataPrey.1, REML=F)



summary(modelPREY1)
sjPlot::tab_df(prey1.ordered[1:10,],
               title = "Physiological Delay Top 10 Models (Prey Availability)", 
               file = "Results/Tables/Prey1Top5.doc")



PREY1 <- list(modelPREY1.1, modelPREY1.2)
PREY1.aic <- round(c(AICc(modelPREY1.1), AICc(modelPREY1.2)), digits=2)
w.PREY1<- prey1.ordered$AICc.Weight[1:2]
PREY1.Best <-  modelPREY1.1

legend.size <- 12
PREY1.plot <-dwplot(PREY1, 
                      vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       allSmolt="Chinook Salmon \n Smolt Production",
                       HakeBiomass="Hake \n Spawning \n Biomass",
                       Chinook = "Chinook Salmon \n Escapement", 
                       Herring.Biomass="Herring \n Spawning \n Biomass"
  )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("C. Physiological Delay") +
  scale_color_manual(name="Model AIC", labels = PREY1.aic, values=hcl.colors(3, palette = "Greens", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY1.plot



###### PREY Lag 2 ####
model.selectionPREY.2 <- model.selection.PREY(dataPrey.2, n, dataPrey.2$TP)
prey2<-data.frame(model.selectionPREY.2)
prey2.ordered <- prey2[order(prey2$AICc),]
prey2.ordered[1:6,]
summary(modelPREY2)
sjPlot::tab_df(prey2.ordered[1:10,],
               title = "1-Year Ecological Top 10 Models (Prey Availability)", 
               file = "Results/Tables/Prey2Top5.doc")

modelPREY2.1<-lmer(TP~Location.2+allSmolt+Chinook+HakeBiomass+(1|AA), data=dataPrey.2, REML=F)
modelPREY2.2<-lmer(TP~Location.2+allSmolt+HakeBiomass+(1|AA), data=dataPrey.2, REML=F)
modelPREY2.3<-lmer(TP~Location.2+allSmolt+Chinook+(1|AA), data=dataPrey.2, REML=F)
modelPREY2.4<-lmer(TP~Location.2+allSmolt+Chinook+HakeBiomass+(1|AA), data=dataPrey.2, REML=F)


PREY2 <- list(modelPREY2.1, modelPREY2.2,
               modelPREY2.3, modelPREY2.4)
PREY2.aic <- round(c(AICc(modelPREY2.1), AICc(modelPREY2.2),
                     AICc(modelPREY2.3), AICc(modelPREY2.4)), digits=2)
w.PREY2<- prey2.ordered$AICc.Weight[1:4]
PREY2.Best <-  modelPREY2.1

legend.size <- 12
PREY2.plot <-dwplot(PREY2, 
                    vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       allSmolt="Chinook Salmon \n Smolt Production",
                       HakeBiomass="Hake \n Spawning \n Biomass",
                       Chinook = "Chinook Salmon \n Escapement", 
                       Herring.Biomass="Herring \n Spawning \n Biomass"
  )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("B. 1-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = PREY2.aic, values=hcl.colors(6, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY2.plot





##### Prey lag 3####
model.selectionPREY.3 <- model.selection.PREY(dataPrey.3, n, dataPrey.3$TP)
prey3<-data.frame(model.selectionPREY.3)
prey3.ordered <- prey3[order(prey3$AICc),]
prey3.ordered[1:7,]
sjPlot::tab_df(prey3.ordered[1:10,],
               title = "2-Year Ecological Top 10 Models (Prey Availability)", 
               file = "Results/Tables/Prey3Top5.doc")

modelPREY3.1<-lmer(TP~Location.2+allSmolt+Chinook+HakeBiomass+(1|AA), data=dataPrey.3, REML=F)
modelPREY3.2<-lmer(TP~Location.2+allSmolt+Chinook+(1|AA), data=dataPrey.3, REML=F)
modelPREY3.3<-lmer(TP~Location.2+HakeBiomass+Chinook+(1|AA), data=dataPrey.3, REML=F)


PREY3 <- list(modelPREY3.1, modelPREY3.2,
              modelPREY3.3)
PREY3.aic <- round(c(AICc(modelPREY3.1), AICc(modelPREY3.2),
                     AICc(modelPREY3.3)), digits=2)
w.PREY3<- prey3.ordered$AICc.Weight[1:3]
PREY3.Best <-  modelPREY3.1

legend.size <- 12
PREY3.plot <-dwplot(PREY3, 
                    vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       allSmolt="Chinook Salmon \n Smolt Production",
                       HakeBiomass="Hake \n Spawning \n Biomass",
                        Chinook = "Chinook Salmon \n Escapement", 
                       Herring.Biomass="Herring \n Spawning \n Biomass"
                       )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("A. 2-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = PREY3.aic, values=hcl.colors(4, palette = "Magenta", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY3.plot


PREY1.plot
PREY2.plot
PREY3.plot
####FULL PLOT####


CLIM.plot<-ggarrange(CLIM3.plot,CLIM2.plot,CLIM1.plot,
                     ncol = 3, nrow = 1,align="hv")+
  theme(plot.margin = margin(0, 0, 0, -.5, "cm"),
        panel.margin=unit(c(0,0,0,0), "cm"))

CLIM.plot<-annotate_figure(CLIM.plot,
                           top = text_grob("Climate Models", color = "black", face = "bold", size = 14)
)

CLIM.plot


pdf(file="Results/Figures/HCoefPlotCLIMATE.pdf", width=12, height=6)
CLIM.plot
dev.off()



PREY.plot<-ggarrange(PREY3.plot,PREY2.plot,PREY1.plot,
                     ncol = 3, nrow = 1,align="hv")+
  theme(plot.margin = margin(0, 0, 0, -.5, "cm"),
        panel.margin=unit(c(0,0,0,0), "cm"))

PREY.plot<-annotate_figure(PREY.plot,
                           top = text_grob("Food Web Models", color = "black", face = "bold", size = 14)
)

PREY.plot


pdf(file="Results/Figures/HCoefPlotFOODWEB.pdf", width=12, height=6)
PREY.plot
dev.off()



############### Residual Time Series #################
palette(c("#ED90A4", "#C0AB52","#28BBD7", "#4FBF85",  "#C699E7"))
palette()

cols<-c('#7EBA68','#CCA65A','#D494E1','#6FB1E7')
palette(c('#7EBA68','#00C1B2','#CCA65A','#6FB1E7','#D494E1'))

dataCLIM.1<- cbind(dataCLIM.1,residuals = residuals(modelCLIM1.1))
dataCLIM.1$AA<- factor(dataCLIM.1$AA, levels = c('ALA','GLU',
                                                 'PRO', 'VAL', 'ASP'),
                       labels = c("Alanine", "Glutamic Acid", "Proline", "Valine*", "Aspartic Acid")
)


clim1_TS <- ggplot(dataCLIM.1, aes(x = Year, y = residuals)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 6), aes(color = AA))+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA~.,labeller = labeller(dataCLIM.1$AA))+
  scale_colour_manual(values = cols)+
  labs(title = "A. Physological Delay Model", y="Residuals", x="Year")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
clim1_TS

dataCLIM.2<- cbind(dataCLIM.2,residuals = residuals(modelCLIM2.1))
dataCLIM.2$AA<- factor(dataCLIM.2$AA, levels = c('ALA','GLU',
                                                 'PRO', 'VAL', 'ASP'),
                       labels = c("Alanine", "Glutamic Acid", "Proline", "Valine*", "Aspartic Acid"))


clim2_TS <- ggplot(dataCLIM.2, aes(x = Year, y = residuals)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 6), aes(color = AA))+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA~.,labeller = labeller(dataCLIM.2$AA))+
  scale_colour_manual(values = cols)+
  labs(title = "B. 1-Year Ecological Delay Model", y="Residuals", x="Year")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
clim2_TS

dataCLIM.3<- cbind(dataCLIM.3,residuals = residuals(modelCLIM3.1))
dataCLIM.3$AA<- factor(dataCLIM.3$AA, levels = c('ALA','GLU',
                                                 'PRO', 'VAL', 'ASP'),
                       labels = c("Alanine", "Glutamic Acid", "Proline", "Valine*", "Aspartic Acid"))


clim3_TS <- ggplot(dataCLIM.3, aes(x = Year, y = residuals)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 6), aes(color = AA))+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA~.,labeller = labeller(dataCLIM.3$AA))+
  scale_colour_manual(values = cols)+
  labs(title = "C. 2-Year Ecological Delay Model", y="Residuals", x="Year")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
clim3_TS


ResidualTSOcean<-annotate_figure(ggarrange(clim1_TS, clim2_TS, clim3_TS, ncol = 3, nrow = 1),top = text_grob("Ocean Condition Models", color = "black", face = "bold", size = 14))


pdf(file="Results/Figures/ResidualTSOcean.pdf", width=11, height=6.5)
ResidualTSOcean
dev.off()






dataPrey.1<- cbind(dataPrey.1,residuals = residuals(modelPREY1.1))
dataPrey.1$AA<- factor(dataPrey.1$AA, levels = c('ALA','GLU',
                                                 'PRO', 'VAL', 'ASP'),
                       labels = c("Alanine", "Glutamic Acid", "Proline", "Valine*", "Aspartic Acid")
)


prey1_TS <- ggplot(dataPrey.1, aes(x = Year, y = residuals)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 6), aes(color = AA))+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA~.,labeller = labeller(dataCLIM.1$AA))+
  scale_colour_manual(values = cols)+
  labs(title = "A. Physological Delay Model", y="Residuals", x="Year")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
prey1_TS

dataPrey.2<- cbind(dataPrey.2,residuals = residuals(modelPREY2.1))
dataPrey.2$AA<- factor(dataPrey.2$AA, levels = c('ALA','GLU',
                                                 'PRO', 'VAL', 'ASP'),
                       labels = c("Alanine", "Glutamic Acid", "Proline", "Valine*", "Aspartic Acid")
)


prey2_TS <- ggplot(dataPrey.2, aes(x = Year, y = residuals)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 6), aes(color = AA))+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA~.,labeller = labeller(dataCLIM.1$AA))+
  scale_colour_manual(values = cols)+
  labs(title = "B. 1-Year Ecological Delay Model", y="Residuals", x="Year")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
prey2_TS

dataPrey.3<- cbind(dataPrey.3,residuals = residuals(modelPREY3.1))
dataPrey.3$AA<- factor(dataPrey.3$AA, levels = c('ALA','GLU',
                                                 'PRO', 'VAL', 'ASP'),
                       labels = c("Alanine", "Glutamic Acid", "Proline", "Valine*", "Aspartic Acid")
)


prey3_TS <- ggplot(dataPrey.3, aes(x = Year, y = residuals)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 6), aes(color = AA))+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA~.,labeller = labeller(dataCLIM.1$AA))+
  scale_colour_manual(values = cols)+
  labs(title = "C. 2-Year Ecological Delay Model", y="Residuals", x="Year")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
prey3_TS


ResidualTSPrey<-annotate_figure(ggarrange(prey1_TS, prey2_TS, prey3_TS, ncol = 3, nrow = 1),top = text_grob("Food Web Models", color = "black", face = "bold", size = 14))


pdf(file="Results/Figures/ResidualTSPrey.pdf", width=11, height=6.5)
ResidualTSPrey
dev.off()




###############Residual Plots ##############
library(HLMdiag)
HLMresid(modelFULL, level=1)
col<-c('#CCA65A','#7EBA68','#00C1B2','#6FB1E7')
palette(c('#7EBA68','#00C1B2','#CCA65A','#6FB1E7','#D494E1'))

pdf(file="Results/Figures/ResidualsDiagnosticsCLIM.pdf", width=7, height=10)
par(mfrow=c(3,2))

y<- dataCLIM.1
plot(predict(modelCLIM1.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="A. Physiological Delay Model", col=dataCLIM.1$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelCLIM1.1), residuals(modelCLIM1.1), ylab = "Residuals", xlab = "Predicted",
     main= "Physiological Delay Model",  pch=19, cex=0.5, col=dataCLIM.1$AA, ylim = c(-2,2))
abline(0, 0)

y<- dataCLIM.2
plot(predict(modelCLIM2.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="B. 1-Year Ecological Delay", col=dataCLIM.2$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelCLIM2.1), residuals(modelCLIM2.1), ylab = "Residuals", xlab = "Predicted",
     main= "1-Year Ecological Delay",  pch=19, cex=0.5, col=dataCLIM.2$AA, ylim = c(-2,2))
abline(0, 0)

y<- dataCLIM.3
plot(predict(modelCLIM3.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="C. 2-Year Ecological Delay", col=dataCLIM.3$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelCLIM3.1), residuals(modelCLIM3.1), ylab = "Residuals", xlab = "Predicted",
     main= "2-Year Ecological Delay",  pch=19, cex=0.5, col=dataCLIM.3$AA, ylim = c(-2,2))
abline(0, 0)



dev.off()


pdf(file="Results/Figures/ResidualsDiagnosticsPREY.pdf", width=7, height=10)
par(mfrow=c(3,2))

y<- dataPrey.1
plot(predict(modelPREY1.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="A. Physiological Delay Model", col=dataPrey.1$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelPREY1.1), residuals(modelPREY1.1), ylab = "Residuals", xlab = "Predicted",
     main= "Physiological Delay Model",  pch=19, cex=0.5, col=dataPrey.1$AA, ylim = c(-2,2))
abline(0, 0)

y<- dataPrey.2
plot(predict(modelPREY2.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="B. 1-Year Ecological Delay", col=dataPrey.2$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelPREY2.1), residuals(modelPREY2.1), ylab = "Residuals", xlab = "Predicted",
     main= "1-Year Ecological Delay",  pch=19, cex=0.5, col=dataPrey.2$AA, ylim = c(-2,2))
abline(0, 0)

y<- dataPrey.3
plot(predict(modelPREY3.1), y[,'TP'], pch=19, cex=0.5, ylab = "Observed", xlab = "Predicted",
     main="C. 2-Year Ecological Delay", col=dataPrey.3$AA)
abline(0, 1)
legend('bottomright', legend = c("Alanine", "Glutamic Acid", "Proline", "Valine"), col = c(1,3:5), cex = 0.8, pch = 16)

plot(predict(modelPREY3.1), residuals(modelPREY3.1), ylab = "Residuals", xlab = "Predicted",
     main= "2-Year Ecological Delay",  pch=19, cex=0.5, col=dataPrey.3$AA, ylim = c(-2,2))
abline(0, 0)



dev.off()
#######TABLES######


##### Model Coefficients #####
confint(modelCLIM1.1, level=0.95)
confint(modelCLIM2.1, level=0.95)
confint(modelCLIM3.1, level=0.95)

confint(modelPREY1.1, level=0.95)
summary(modelPREY1.1)

confint(modelPREY2.1, level=0.95)
summary(modelPREY2.1)
