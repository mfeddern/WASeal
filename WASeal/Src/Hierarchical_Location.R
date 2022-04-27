
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
                      #AICc(lmer(y~(1|AA), data=dataframe, REML=F)), #2
                      AICc(lmer(y~PDO+(1|AA), data=dataframe, REML=F)),#3
                      AICc(lmer(y~NPGO+(1|AA), data=dataframe, REML=F)),#4
                      AICc(lmer(y~MEI+(1|AA), data=dataframe, REML=F)),#5
                      AICc(lmer(y~UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#6
                      AICc(lmer(y~PDO+NPGO+(1|AA), data=dataframe, REML=F)),#7
                      AICc(lmer(y~PDO+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#8
                      AICc(lmer(y~UpInAn.45.Spring+NPGO+(1|AA), data=dataframe, REML=F)),#9
                      AICc(lmer(y~MEI+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#10
                      AICc(lmer(y~(1|AA)+WA.SST.Su, data=dataframe, REML=F)), #11
                      AICc(lmer(y~WA.SST.Su+PDO+(1|AA), data=dataframe, REML=F)),#12
                      AICc(lmer(y~WA.SST.Su+NPGO+(1|AA), data=dataframe, REML=F)),#13
                      AICc(lmer(y~WA.SST.Su+MEI+(1|AA), data=dataframe, REML=F)),#14
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#15
                      AICc(lmer(y~WA.SST.Su+PDO+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#16
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Spring+NPGO+(1|AA), data=dataframe, REML=F)),#17
                      AICc(lmer(y~WA.SST.Su+MEI+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#18
                      AICc(lmer(y~(1|AA)+UpInAn.45.Summer, data=dataframe, REML=F)), #19
                      AICc(lmer(y~UpInAn.45.Summer+PDO+(1|AA), data=dataframe, REML=F)),#20
                      AICc(lmer(y~UpInAn.45.Summer+NPGO+(1|AA), data=dataframe, REML=F)),#21
                      AICc(lmer(y~UpInAn.45.Summer+MEI+(1|AA), data=dataframe, REML=F)),#22
                      AICc(lmer(y~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+(1|AA), data=dataframe, REML=F)),#23
                      AICc(lmer(y~(1|AA)+Col.Dis.high, data=dataframe, REML=F)), #24
                      AICc(lmer(y~Col.Dis.high+PDO+(1|AA), data=dataframe, REML=F)),#25
                      AICc(lmer(y~Col.Dis.high+NPGO+(1|AA), data=dataframe, REML=F)),#26
                      AICc(lmer(y~Col.Dis.high+MEI+(1|AA), data=dataframe, REML=F)),#27
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#28
                      AICc(lmer(y~Col.Dis.high+PDO+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#29
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Spring+NPGO+(1|AA), data=dataframe, REML=F)),#30
                      AICc(lmer(y~Col.Dis.high+MEI+UpInAn.45.Spring+(1|AA), data=dataframe, REML=F)),#31
                      AICc(lmer(y~Col.Dis.high+WA.SST.Su+(1|AA), data=dataframe, REML=F)),#32
                      AICc(lmer(y~Col.Dis.high+UpInAn.45.Summer+(1|AA), data=dataframe, REML=F)),#33
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+(1|AA), data=dataframe, REML=F)),#34
                      AICc(lmer(y~WA.SST.Su+UpInAn.45.Summer+Col.Dis.high+(1|AA), data=dataframe, REML=F))#35

                      
                      
                      
  )
  
  model.names <- c("1. Null",
                   #"2. Location Only", 
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

dataCLIM.1.coastal<- subset(dataCLIM.1, Location.2=="Coastal")
model.selectionCLIM.1 <- ModelSelection.CLIM(dataCLIM.1.coastal, n, dataCLIM.1.coastal$TP)
clim.1 <-data.frame(model.selectionCLIM.1)
clim.1.ordered <- clim.1[order(clim.1$AICc),]
clim.1.ordered[1:5,]

modelCLIM1.1.coastal<-lmer(TP~UpInAn.45.Spring+WA.SST.Su+PDO+(1|AA), data=dataCLIM.1.coastal, REML=F)
modelCLIM1.2.coastal<-lmer(TP~UpInAn.45.Summer+NPGO+(1|AA), data=dataCLIM.1.coastal, REML=F)
modelCLIM1.3.coastal<-lmer(TP~Col.Dis.high+MEI+(1|AA), data=dataCLIM.1.coastal, REML=F)
modelCLIM1.4.coastal<-lmer(TP~WA.SST.Su+MEI+(1|AA), data=dataCLIM.1.coastal, REML=F)
modelCLIM1.5.coastal<-lmer(TP~UpInAn.45.Summer+UpInAn.45.Spring+NPGO+(1|AA), data=dataCLIM.1.coastal, REML=F)



dataCLIM.1.inland<- subset(dataCLIM.1, Location.2=="Inland")
model.selectionCLIM.1 <- ModelSelection.CLIM(dataCLIM.1.inland, n, dataCLIM.1.inland$TP)
clim.1 <-data.frame(model.selectionCLIM.1)
clim.1.ordered <- clim.1[order(clim.1$AICc),]
clim.1.ordered[1:5,]

modelCLIM1.1.inland<-lmer(TP~UpInAn.45.Spring+UpInAn.45.Summer+NPGO+(1|AA), data=dataCLIM.1.inland, REML=F)
modelCLIM1.2.inland<-lmer(TP~UpInAn.45.Summer+Col.Dis.high+(1|AA), data=dataCLIM.1.inland, REML=F)
modelCLIM1.3.inland<-lmer(TP~UpInAn.45.Summer+Col.Dis.high+WA.SST.Su+(1|AA), data=dataCLIM.1.inland, REML=F)
modelCLIM1.4.inland<-lmer(TP~UpInAn.45.Summer+NPGO+(1|AA), data=dataCLIM.1.inland, REML=F)
modelCLIM1.5.inland<-lmer(TP~UpInAn.45.Summer+(1|AA), data=dataCLIM.1.inland, REML=F)

#### Plots Climate 1 ####

CLIM1 <- list(modelCLIM1.1.coastal, modelCLIM1.2.coastal,
              modelCLIM1.3.coastal,modelCLIM1.4.coastal,
              modelCLIM1.5.coastal)
CLIM1.aic <- round(c(AICc(modelCLIM1.1.coastal), AICc(modelCLIM1.2.coastal),
                     AICc(modelCLIM1.3.coastal), AICc(modelCLIM1.4.coastal),
                     AICc(modelCLIM1.5.coastal)), digits=2)
w.CLIM1<- clim.1.ordered$AICc.Weight[1:5]
CLIM1.Best <-  modelCLIM1.1.coastal

legend.size <- 12
CLIM1.plot.coastal <-dwplot(CLIM1, 
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

CLIM1.plot.coastal



CLIM1 <- list(modelCLIM1.1.inland, modelCLIM1.2.inland, modelCLIM1.3.inland,
              modelCLIM1.4.inland,modelCLIM1.5.inland)

CLIM1.aic <- round(c(AICc(modelCLIM1.1.inland), AICc(modelCLIM1.2.inland),
                     AICc(modelCLIM1.3.inland), AICc(modelCLIM1.4.inland),
                     AICc(modelCLIM1.5.inland)), digits=2)
w.CLIM1<- clim.1.ordered$AICc.Weight[1:5]
CLIM1.Best <-  modelCLIM1.1.inland
legend.size <- 12
CLIM1.plot.inland <-dwplot(CLIM1, 
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

CLIM1.plot.inland


##### CLIMATE 2#####

dataCLIM.2.coastal<- subset(dataCLIM.2, Location.2=="Coastal")
model.selectionCLIM.2 <- ModelSelection.CLIM(dataCLIM.2.coastal, n, dataCLIM.2.coastal$TP)
clim.2 <-data.frame(model.selectionCLIM.2)
clim.2.ordered <- clim.2[order(clim.2$AICc),]
clim.2.ordered[1:5,]


modelCLIM2.1.coastal<-lmer(TP~UpInAn.45.Summer+WA.SST.Su+Col.Dis.high+(1|AA), data=dataCLIM.2.coastal, REML=F)
modelCLIM2.2.coastal<-lmer(TP~UpInAn.45.Summer+WA.SST.Su+(1|AA), data=dataCLIM.2.coastal, REML=F)
modelCLIM2.3.coastal<-lmer(TP~NPGO+WA.SST.Su+(1|AA), data=dataCLIM.2.coastal, REML=F)
modelCLIM2.4.coastal<-lmer(TP~Col.Dis.high+WA.SST.Su+(1|AA), data=dataCLIM.2.coastal, REML=F)
modelCLIM2.5.coastal<-lmer(TP~NPGO+(1|AA), data=dataCLIM.2.coastal, REML=F)


dataCLIM.2.inland<- subset(dataCLIM.2, Location.2=="Inland")
model.selectionCLIM.2 <- ModelSelection.CLIM(dataCLIM.2.inland, n, dataCLIM.2.inland$TP)
clim.2 <-data.frame(model.selectionCLIM.2)
clim.2.ordered <- clim.2[order(clim.2$AICc),]
clim.2.ordered[1:5,]

modelCLIM2.1.inland<-lmer(TP~UpInAn.45.Spring+WA.SST.Su+NPGO+(1|AA), data=dataCLIM.2.inland, REML=F)
modelCLIM2.2.inland<-lmer(TP~UpInAn.45.Spring+WA.SST.Su+(1|AA), data=dataCLIM.2.inland, REML=F)
modelCLIM2.3.inland<-lmer(TP~UpInAn.45.Spring+PDO+WA.SST.Su+(1|AA), data=dataCLIM.2.inland, REML=F)
modelCLIM2.4.inland<-lmer(TP~UpInAn.45.Spring+MEI+WA.SST.Su+(1|AA), data=dataCLIM.2.inland, REML=F)
modelCLIM2.5.inland<-lmer(TP~NPGO+WA.SST.Su+(1|AA), data=dataCLIM.2.inland, REML=F)

#### Plots Climate 2 ####

CLIM2 <- list(modelCLIM2.1.coastal, modelCLIM2.2.coastal,
              modelCLIM2.3.coastal,modelCLIM2.4.coastal,modelCLIM2.5.coastal)
CLIM2.aic <- round(c(AICc(modelCLIM2.1.coastal), AICc(modelCLIM2.2.coastal),
                     AICc(modelCLIM2.3.coastal), AICc(modelCLIM2.4.coastal),
                     AICc(modelCLIM2.5.coastal)), digits=2)
w.CLIM2<- clim.2.ordered$AICc.Weight[1:5]
CLIM2.Best <-  modelCLIM2.1.coastal


legend.size <- 12
CLIM2.plot.coastal <-dwplot(CLIM2, 
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
  ggtitle(aes(label="B. 1-Year Ecological Delay", fontface=plain)) +
  scale_color_manual(name="Model AIC", labels = CLIM1.aic, values=hcl.colors(7, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

CLIM2.plot.coastal



CLIM2 <- list(modelCLIM2.1.inland, modelCLIM2.2.inland,modelCLIM2.3.inland,
              modelCLIM2.4.inland,modelCLIM2.5.inland)
CLIM2.aic <- round(c(AICc(modelCLIM2.1.inland), AICc(modelCLIM2.2.inland),
                     AICc(modelCLIM2.3.inland), AICc(modelCLIM2.4.inland),
                     AICc(modelCLIM2.5.inland)), digits=2)
w.CLIM2<- clim.2.ordered$AICc.Weight[1:5]
CLIM2.Best <-  modelCLIM2.1.inland

legend.size <- 12
CLIM2.plot.inland <-dwplot(CLIM2, 
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
  ggtitle(aes(label="B. 1-Year Ecological Delay", fontface=plain)) +
  scale_color_manual(name="Model AIC", labels = CLIM1.aic, values=hcl.colors(7, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

CLIM2.plot.inland


#### CLIMATE 3 #####

dataCLIM.3.coastal<- subset(dataCLIM.3, Location.2=="Coastal")
model.selectionCLIM.3 <- ModelSelection.CLIM(dataCLIM.3.coastal, n, dataCLIM.3.coastal$TP)
clim.3 <-data.frame(model.selectionCLIM.3)
clim.3.ordered <- clim.3[order(clim.3$AICc),]
clim.3.ordered[1:5,]

modelCLIM3.1.coastal<-lmer(TP~UpInAn.45.Summer+Col.Dis.high+(1|AA), data=dataCLIM.3.coastal, REML=F)
modelCLIM3.2.coastal<-lmer(TP~UpInAn.45.Summer+WA.SST.Su+Col.Dis.high+(1|AA), data=dataCLIM.3.coastal, REML=F)
modelCLIM3.3.coastal<-lmer(TP~WA.SST.Su+Col.Dis.high+(1|AA), data=dataCLIM.3.coastal, REML=F)
modelCLIM3.4.coastal<-lmer(TP~Col.Dis.high+UpInAn.45.Spring+PDO+(1|AA), data=dataCLIM.3.coastal, REML=F)
modelCLIM3.5.coastal<-lmer(TP~Col.Dis.high+(1|AA), data=dataCLIM.3.coastal, REML=F)


dataCLIM.3.inland<- subset(dataCLIM.3, Location.2=="Inland")
model.selectionCLIM.3 <- ModelSelection.CLIM(dataCLIM.3.inland, n, dataCLIM.3.inland$TP)
clim.3 <-data.frame(model.selectionCLIM.3)
clim.3.ordered <- clim.3[order(clim.3$AICc),]
clim.3.ordered[1:5,]

modelCLIM3.1.inland<-lmer(TP~UpInAn.45.Spring+Col.Dis.high+(1|AA), data=dataCLIM.3.inland, REML=F)
modelCLIM3.2.inland<-lmer(TP~WA.SST.Su+Col.Dis.high+(1|AA), data=dataCLIM.3.inland, REML=F)
modelCLIM3.3.inland<-lmer(TP~UpInAn.45.Spring+(1|AA), data=dataCLIM.3.inland, REML=F)
modelCLIM3.4.inland<-lmer(TP~Col.Dis.high+NPGO+(1|AA), data=dataCLIM.3.inland, REML=F)
modelCLIM3.5.inland<-lmer(TP~Col.Dis.high+(1|AA), data=dataCLIM.3.inland, REML=F)

#### Plots Climate 3 ####

CLIM3 <- list(modelCLIM3.1.coastal, modelCLIM3.2.coastal,modelCLIM3.3.coastal,
              modelCLIM3.4.coastal,modelCLIM3.5.coastal)
CLIM3.aic <- round(c(AICc(modelCLIM3.1.coastal), AICc(modelCLIM3.2.coastal),
                     AICc(modelCLIM3.3.coastal), AICc(modelCLIM3.4.coastal),
                     AICc(modelCLIM3.5.coastal)), digits=2)
w.CLIM3<- clim.3.ordered$AICc.Weight[1:5]
CLIM3.Best <-  modelCLIM3.1.coastal

legend.size <- 12
CLIM3.plot.coastal <-dwplot(CLIM3, 
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
  ggtitle(aes(label="A. 2-Year Ecological Delay", fontface=plain)) +
  scale_color_manual(name="Model AIC", labels = CLIM1.aic, values=hcl.colors(7, palette = "Magenta", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

CLIM3.plot.coastal

CLIM3 <- list(modelCLIM3.1.inland, modelCLIM3.2.inland, modelCLIM3.3.inland, modelCLIM3.4.inland, modelCLIM3.5.inland)
CLIM3.aic <- round(c(AICc(modelCLIM3.1.inland), AICc(modelCLIM3.2.inland), AICc(modelCLIM3.3.inland),
                     AICc(modelCLIM3.4.inland),AICc(modelCLIM3.5.inland)), digits=2)
w.CLIM3<- clim.3.ordered$AICc.Weight[1:5]
CLIM3.Best <-  modelCLIM3.1.inland

legend.size <- 12
CLIM3.plot.inland <-dwplot(CLIM3, 
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
  ggtitle(aes(label="A. 2-Year Ecological Delay", fontface=plain)) +
  scale_color_manual(name="Model AIC", labels = CLIM1.aic, values=hcl.colors(7, palette = "Magenta", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

CLIM3.plot.inland

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
                      #AICc(lmer(y~(1|AA), data=dataframe, REML=F)), #2 Location Only
                      AICc(lmer(y~Herring.Biomass+(1|AA), data=dataframe, REML=F)), #3. Herring Biomass
                      AICc(lmer(y~Chinook+(1|AA), data=dataframe, REML=F)),#4. Chinook Escapements
                      AICc(lmer(y~allSmolt+(1|AA), data=dataframe, REML=F)), #5. All smolts
                      AICc(lmer(y~HakeBiomass+(1|AA), data=dataframe, REML=F)), #6. Hake Biomass
                      AICc(lmer(y~Herring.Biomass+Chinook+(1|AA), data=dataframe, REML=F)),#7. Herring Biomass, Chinook Escapements
                      AICc(lmer(y~Herring.Biomass+HakeBiomass+(1|AA), data=dataframe, REML=F)),#8. Herring Biomass, Hake Biomass
                      AICc(lmer(y~Herring.Biomass+allSmolt+(1|AA), data=dataframe, REML=F)),#9. Herring Biomass, all smolt
                      AICc(lmer(y~HakeBiomass+Chinook+(1|AA), data=dataframe, REML=F)),#10. Hake Bimass, Chinook Escapement
                      AICc(lmer(y~allSmolt+Chinook+(1|AA), data=dataframe, REML=F)),#11. all smolt, Chinook Escapements
                      AICc(lmer(y~allSmolt+HakeBiomass+(1|AA), data=dataframe, REML=F)),#12. all smolt, Hake Biomass
                      AICc(lmer(y~allSmolt+Chinook+HakeBiomass+(1|AA), data=dataframe, REML=F)),#13. all smolt, Chinook escapent, Hake Biomas
                      AICc(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+(1|AA), data=dataframe, REML=F)),#14. Herring Biomass, all smolts, Hake Biomass
                      AICc(lmer(y~Herring.Biomass+allSmolt+Chinook+(1|AA), data=dataframe, REML=F)),#15.Herring biomass, all smolts, Chinook escapements
                      AICc(lmer(y~Herring.Biomass+Chinook+HakeBiomass+(1|AA), data=dataframe, REML=F)), #16. Herring Biomass, Chinook escapements, Hake biomass
                      AICc(lmer(y~HarborSeal+(1|AA), data=dataframe, REML=F)), #17. Harbor seal 
                      AICc(lmer(y~HarborSeal+Herring.Biomass+(1|AA), data=dataframe, REML=F)),#18. Harbor seal, herring biomass
                      AICc(lmer(y~HarborSeal+Chinook+(1|AA), data=dataframe, REML=F)), #19. Harbor seal, Chinook escapement
                      AICc(lmer(y~HarborSeal+allSmolt+(1|AA), data=dataframe, REML=F)),#20. Harbor seal, all smolt
                      AICc(lmer(y~HarborSeal+HakeBiomass+(1|AA), data=dataframe, REML=F)),#21. Harbor seal, hake biomass
                      AICc(lmer(y~HarborSeal+Herring.Biomass+Chinook+(1|AA), data=dataframe, REML=F)),#22. Harbor seal, herring biomass, Chinook escapement
                      AICc(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+(1|AA), data=dataframe, REML=F)),#23. Harbor seal, herring biomass, hakr biomass
                      AICc(lmer(y~HarborSeal+Herring.Biomass+allSmolt+(1|AA), data=dataframe, REML=F)),#24. Harbor seal, herring biomass, all smolts
                      AICc(lmer(y~HarborSeal+HakeBiomass+Chinook+(1|AA), data=dataframe, REML=F)),#25. Harbor seal, hake biomass, Chinool escapement
                      AICc(lmer(y~HarborSeal+allSmolt+HakeBiomass+(1|AA), data=dataframe, REML=F))#26. Harbor seal, all smolts, hake biomass
                    #  AICc(lmer(y~Herring.Biomass*HakeBiomass+Location.2+(1|AA), data=dataframe, REML=F)) #27. Hake:Herring biomass
                      
                      
                      
  )
  
  model.names <- c("1. Null",
                  # "2. Location Only",
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

dataPREY.1.coastal<- subset(dataPrey.1, Location.2=="Coastal")
model.selectionPREY.1 <- model.selection.PREY(dataPREY.1.coastal, n, dataPREY.1.coastal$TP)
prey.1 <-data.frame(model.selectionPREY.1)
prey.1.ordered <- prey.1[order(prey.1$AICc),]
prey.1.ordered[1:10,]
modelPREY1.1.coastal<-lmer(TP~HarborSeal+Herring.Biomass+HakeBiomass+(1|AA), data=dataPREY.1.coastal, REML=F)
modelPREY1.2.coastal<-lmer(TP~allSmolt+HakeBiomass+(1|AA), data=dataPREY.1.coastal, REML=F)
modelPREY1.3.coastal<-lmer(TP~allSmolt+Herring.Biomass+HakeBiomass+(1|AA), data=dataPREY.1.coastal, REML=F)
modelPREY1.4.coastal<-lmer(TP~allSmolt+(1|AA), data=dataPREY.1.coastal, REML=F)
modelPREY1.5.coastal<-lmer(TP~HarborSeal+allSmolt+HakeBiomass+(1|AA), data=dataPREY.1.coastal, REML=F)

dataPREY.1.inland<- subset(dataPrey.1, Location.2=="Inland")
model.selectionPREY.1 <- model.selection.PREY(dataPREY.1.inland, n, dataPREY.1.inland$TP)
prey.1 <-data.frame(model.selectionPREY.1)
prey.1.ordered <- prey.1[order(prey.1$AICc),]
prey.1.ordered[1:10,]

modelPREY1.1.inland<-lmer(TP~HakeBiomass+(1|AA), data=dataPREY.1.inland, REML=F)
modelPREY1.2.inland<-lmer(TP~allSmolt+HakeBiomass+(1|AA), data=dataPREY.1.inland, REML=F)
modelPREY1.3.inland<-lmer(TP~Chinook+HakeBiomass+(1|AA), data=dataPREY.1.inland, REML=F)
modelPREY1.4.inland<-lmer(TP~HarborSeal+HakeBiomass+(1|AA), data=dataPREY.1.inland, REML=F)
modelPREY1.5.inland<-lmer(TP~Herring.Biomass+HakeBiomass+(1|AA), data=dataPREY.1.inland, REML=F)


#### Plots Prey 1 ####
PREY1.coastal <- list(modelPREY1.1.coastal, modelPREY1.2.coastal, modelPREY1.3.coastal, modelPREY1.4.coastal, modelPREY1.5.coastal)
PREY1.aic <- round(c(AICc(modelPREY1.1.coastal), AICc(modelPREY1.2.coastal),
                     AICc(modelPREY1.3.coastal), AICc(modelPREY1.4.coastal),
                     AICc(modelPREY1.5.coastal)), digits=2)
w.PREY1<- prey1.ordered$AICc.Weight[1:2]
PREY1.Best <-  modelPREY1.1

legend.size <- 12
PREY1.plot.coastal <-dwplot(PREY1.coastal, 
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
  scale_color_manual(name="Model AIC", labels = PREY1.aic, values=hcl.colors(5, palette = "Greens", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY1.plot.coastal

PREY1.inland <- list(modelPREY1.1.inland, modelPREY1.2.inland, modelPREY1.3.inland, modelPREY1.4.inland, modelPREY1.5.inland)
PREY1.aic <- round(c(AICc(modelPREY1.1.inland), AICc(modelPREY1.2.inland),
                     AICc(modelPREY1.3.inland), AICc(modelPREY1.4.inland),
                     AICc(modelPREY1.5.inland)), digits=2)
w.PREY1<- prey1.ordered$AICc.Weight[1:2]
PREY1.Best <-  modelPREY1.1

legend.size <- 12
PREY1.plot.inland <-dwplot(PREY1.inland, 
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
  scale_color_manual(name="Model AIC", labels = PREY1.aic, values=hcl.colors(5, palette = "Greens", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY1.plot.inland

#######PREY LAG 2 #######

dataPREY.2.coastal<- subset(dataPrey.2, Location.2=="Coastal")
model.selectionPREY.2 <- model.selection.PREY(dataPREY.2.coastal, n, dataPREY.2.coastal$TP)
prey.2 <-data.frame(model.selectionPREY.2)
prey.2.ordered <- prey.2[order(prey.2$AICc),]
prey.2.ordered[1:5,]

modelPREY2.1.coastal<-lmer(TP~Chinook+(1|AA), data=dataPREY.2.coastal, REML=F)
modelPREY2.2.coastal<-lmer(TP~Chinook+HakeBiomass+(1|AA), data=dataPREY.2.coastal, REML=F)
modelPREY2.3.coastal<-lmer(TP~HarborSeal+Herring.Biomass+HakeBiomass+(1|AA), data=dataPREY.2.coastal, REML=F)
modelPREY2.4.coastal<-lmer(TP~HarborSeal+Chinook+HakeBiomass+(1|AA), data=dataPREY.2.coastal, REML=F)
modelPREY2.5.coastal<-lmer(TP~HarborSeal+Chinook+(1|AA), data=dataPREY.2.coastal, REML=F)


dataPREY.2.inland<- subset(dataPrey.2, Location.2=="Inland")
model.selectionPREY.2 <- model.selection.PREY(dataPREY.2.inland, n, dataPREY.2.inland$TP)
prey.2 <-data.frame(model.selectionPREY.2)
prey.2.ordered <- prey.2[order(prey.2$AICc),]
prey.2.ordered[1:5,]

modelPREY2.1.inland<-lmer(TP~allSmolt+(1|AA), data=dataPREY.2.inland, REML=F)
modelPREY2.2.inland<-lmer(TP~allSmolt+HarborSeal+(1|AA), data=dataPREY.2.inland, REML=F)
modelPREY2.3.inland<-lmer(TP~Chinook+allSmolt+(1|AA), data=dataPREY.2.inland, REML=F)
modelPREY2.4.inland<-lmer(TP~allSmolt+HakeBiomass+(1|AA), data=dataPREY.2.inland, REML=F)
modelPREY2.5.inland<-lmer(TP~Herring.Biomass+allSmolt+(1|AA), data=dataPREY.2.inland, REML=F)

#### Plots Prey 2 ####
PREY2.coastal <- list(modelPREY2.1.coastal, modelPREY2.2.coastal, modelPREY2.3.coastal, modelPREY2.4.coastal, modelPREY2.5.coastal)
PREY2.aic <- round(c(AICc(modelPREY2.1.coastal), AICc(modelPREY2.2.coastal),
                     AICc(modelPREY2.3.coastal), AICc(modelPREY2.4.coastal),
                     AICc(modelPREY2.5.coastal)), digits=2)
w.PREY2<- prey2.ordered$AICc.Weight[1:2]
PREY2.Best <-  modelPREY2.1

legend.size <- 12
PREY2.plot.coastal <-dwplot(PREY2.coastal, 
                            vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       allSmolt="Chinook Salmon \n Smolt Production",
                       HakeBiomass="Hake \n Spawning \n Biomass",
                       Chinook = "Chinook Salmon \n Escapement", 
                       Herring.Biomass="Herring \n Spawning \n Biomass"
  )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("E. 1-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = PREY1.aic, values=hcl.colors(5, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY2.plot.coastal

PREY2.inland <- list(modelPREY2.1.inland, modelPREY2.2.inland, modelPREY2.3.inland, modelPREY2.4.inland, modelPREY2.5.inland)
PREY2.aic <- round(c(AICc(modelPREY2.1.inland), AICc(modelPREY2.2.inland),
                     AICc(modelPREY2.3.inland), AICc(modelPREY2.4.inland),
                     AICc(modelPREY2.5.inland)), digits=2)
w.PREY2<- prey2.ordered$AICc.Weight[1:2]
PREY2.Best <-  modelPREY2.1

legend.size <- 12
PREY2.plot.inland <-dwplot(PREY2.inland, 
                           vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       allSmolt="Chinook Salmon \n Smolt Production",
                       HakeBiomass="Hake \n Spawning \n Biomass",
                       Chinook = "Chinook Salmon \n Escapement", 
                       Herring.Biomass="Herring \n Spawning \n Biomass"
  )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("E. 1-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = PREY1.aic, values=hcl.colors(5, palette = "Teal", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY2.plot.inland


#######PREY LAG 3 #######

dataPREY.3.coastal<- subset(dataPrey.3, Location.2=="Coastal")
model.selectionPREY.3 <- model.selection.PREY(dataPREY.3.coastal, n, dataPREY.3.coastal$TP)
prey.3 <-data.frame(model.selectionPREY.3)
prey.3.ordered <- prey.3[order(prey.3$AICc),]
prey.3.ordered[1:5,]

modelPREY3.1.coastal<-lmer(TP~Herring.Biomass+(1|AA), data=dataPREY.3.coastal, REML=F)
modelPREY3.2.coastal<-lmer(TP~HarborSeal+Herring.Biomass+(1|AA), data=dataPREY.3.coastal, REML=F)
modelPREY3.3.coastal<-lmer(TP~Herring.Biomass+HakeBiomass+(1|AA), data=dataPREY.3.coastal, REML=F)
modelPREY3.4.coastal<-lmer(TP~allSmolt+Chinook+HakeBiomass+(1|AA), data=dataPREY.3.coastal, REML=F)
modelPREY3.5.coastal<-lmer(TP~Herring.Biomass+Chinook+(1|AA), data=dataPREY.3.coastal, REML=F)


dataPREY.3.inland<- subset(dataPrey.3, Location.2=="Inland")
model.selectionPREY.3 <- model.selection.PREY(dataPREY.3.inland, n, dataPREY.3.inland$TP)
prey.3 <-data.frame(model.selectionPREY.3)
prey.3.ordered <- prey.3[order(prey.3$AICc),]
prey.3.ordered[1:10,]

modelPREY3.1.inland<-lmer(TP~allSmolt+(1|AA), data=dataPREY.3.inland, REML=F)
modelPREY3.2.inland<-lmer(TP~allSmolt+Chinook+(1|AA), data=dataPREY.3.inland, REML=F)
modelPREY3.3.inland<-lmer(TP~Herring.Biomass+HakeBiomass+allSmolt+(1|AA), data=dataPREY.3.inland, REML=F)
modelPREY3.4.inland<-lmer(TP~allSmolt+Herring.Biomass+(1|AA), data=dataPREY.3.inland, REML=F)
modelPREY3.5.inland<-lmer(TP~HarborSeal+Chinook+HakeBiomass+(1|AA), data=dataPREY.3.inland, REML=F)


#### Plots Prey 3 ####

PREY3.coastal <- list(modelPREY3.1.coastal, modelPREY3.2.coastal, modelPREY3.3.coastal, modelPREY3.4.coastal, modelPREY3.5.coastal)
PREY3.aic <- round(c(AICc(modelPREY3.1.coastal), AICc(modelPREY3.2.coastal),
                     AICc(modelPREY3.3.coastal), AICc(modelPREY3.4.coastal),
                     AICc(modelPREY3.5.coastal)), digits=2)
w.PREY3<- prey3.ordered$AICc.Weight[1:2]
PREY3.Best <-  modelPREY3.1

legend.size <- 12
PREY3.plot.coastal <-dwplot(PREY3.coastal, 
                            vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       allSmolt="Chinook Salmon \n Smolt Production",
                       HakeBiomass="Hake \n Spawning \n Biomass",
                       Chinook = "Chinook Salmon \n Escapement", 
                       Herring.Biomass="Herring \n Spawning \n Biomass"
  )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("F. 2-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = PREY1.aic, values=hcl.colors(5, palette = "Magenta", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY3.plot.coastal

PREY3.inland <- list(modelPREY3.1.inland, modelPREY3.2.inland, modelPREY3.3.inland, modelPREY3.4.inland, modelPREY3.5.inland)
PREY3.aic <- round(c(AICc(modelPREY3.1.inland), AICc(modelPREY3.2.inland),
                     AICc(modelPREY3.3.inland), AICc(modelPREY3.4.inland),
                     AICc(modelPREY3.5.inland)), digits=2)
w.PREY3<- prey3.ordered$AICc.Weight[1:2]
PREY3.Best <-  modelPREY3.1


PREY3.plot.inland <-dwplot(PREY3.inland, 
                           vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) %>% # plot line at zero _behind_ coefs
  relabel_predictors(c(Location.2Inland = "Subregion \n (Salish Sea)",
                       allSmolt="Chinook Salmon \n Smolt Production",
                       HakeBiomass="Hake \n Spawning \n Biomass",
                       Chinook = "Chinook Salmon \n Escapement", 
                       Herring.Biomass="Herring \n Spawning \n Biomass"
  )) +
  theme_bw() + xlab("") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("F. 2-Year Ecological Delay") +
  scale_color_manual(name="Model AIC", labels = PREY1.aic, values=hcl.colors(5, palette = "Magenta", alpha = NULL, rev = FALSE, fixup = TRUE))+
  theme(plot.title = element_text(face="bold", hjust=0.5, size=legend.size+2),
        plot.margin = margin(0.25, 0, 0, 0, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title.align = .5,
        legend.text=element_text(size=legend.size),
        axis.text.x=element_text(size=legend.size),
        axis.text.y=element_text(size=legend.size),
        plot.subtitle = element_text(hjust = 0.5, size=13))

PREY3.plot.inland


####FULL PLOT####


FULL.plot<-ggarrange(CLIM3.plot.coastal,CLIM2.plot.coastal,CLIM1.plot.coastal,
                    PREY3.plot.coastal,PREY2.plot.coastal,PREY1.plot.coastal,
                     ncol = 3, nrow = 2,align="hv")+
  theme(plot.margin = margin(0, 0, 0, -.5, "cm"),
        panel.margin=unit(c(0,0,0,0), "cm"))

FULL.plot<-annotate_figure(FULL.plot,
                           bottom = text_grob("Coefficient Estimate", color = "black", face = "plain", size = 14)
)

FULL.plot


pdf(file="Results/Figures/HCoefPlotFULL.coastal.pdf", width=12, height=9)
FULL.plot
dev.off()



FULL.plot<-ggarrange(CLIM3.plot.inland,CLIM2.plot.inland,CLIM1.plot.inland,
                     PREY3.plot.inland,PREY2.plot.inland,PREY1.plot.inland,
                     ncol = 3, nrow = 2,align="hv")+
  theme(plot.margin = margin(0, 0, 0, -.5, "cm"),
        panel.margin=unit(c(0,0,0,0), "cm"))

FULL.plot<-annotate_figure(FULL.plot,
                           bottom = text_grob("Coefficient Estimate", color = "black", face = "plain", size = 14)
)

FULL.plot


pdf(file="Results/Figures/HCoefPlotFULL.inland.pdf", width=12, height=9)
FULL.plot
dev.off()
