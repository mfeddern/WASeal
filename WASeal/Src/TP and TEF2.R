rm(list = ls())
library(mgcViz)
library(MARSS)
library(forecast)
library(EBImage)
library(dplyr)
library(lme4)
library(tidyr)
library(dplyr)
#install.packages("digitize")
#install.packages("mgcv")
library(mgcv)
#data <- read.csv("Data/Compiled/WASealAAandTP2.csv")
data <- read.csv("Data/Compiled/SI Data.csv")
#install.packages("visreg")
library(visreg)
#data <- read.csv("SortedTP.csv")
#dev.off()

library(tidyverse)
library(caret)
library(glmnet)
library(Momocs)
library(dplyr)
library(glinternet)
# Loading the data
#data(swiss)

########################     Calculating TEFs                      ############################
########################     Summary Statistics                      ############################


par(mfrow=c(1,1), oma=c(0,0,0,0))
plot(data$years, data$TPmean, pch=16, ylab="Trophic Position", xlab="Year",ylim=c(2,5.5))
summary(lm(TP.GLU~years, data=data))
summary(gam(TP.GLU~s(years)+Sex, data=data)) #gam with a smooth term to years, two level factors with a plus after years
#lets relationship be non linear, 
summary(gam(TP.GLU~s(years, by=Sex), data=data)) #gam with a smooth term to years, two level factors with a plus after years

Location <-subset(data, Location.2=="Inland"|Location.2=="Coastal")
In <- subset(data, Location.2=="Inland")
Co <- subset(data, Location.2=="Coastal")

Sex<-subset(data, Sex=="M"|Sex=="F")
M <- subset(data, Sex=="M")
FM <- subset(data, Sex=="F")

summary(gam(TP.GLU~s(years)+Sex, data=Sex)) #gam with a smooth term to year, two level factors with a plus after year
#lets relationship be non linear, 
summary(gam(TP.GLU~s(years, by=Sex), data=Sex)) #gam with a smooth term to year, two level factors with a plus after year



plot(In$years, In$TPmean, pch=16, col='blue',ylab="Trophic Position",xlab="Year", ylim=c(2,5.5))
points(Co$years, Co$TPmean, pch=16, col='green3')
legend(1925, 5.63, legend=c("Salish Sea", "Coastal"),
       col=c("blue", "green3"), pch=16, cex=0.8)
summary(lm(TPmean~years+Location.2, data=Location))
t.test(TPmean~Location.2, data=Location) #p=0.01252




summary(gam(TP.GLU~s(Year)+Location.2, data=Location)) #gam with a smooth term to year, two level factors with a plus after year
AIC(gam(TP.GLU~s(Year)+Location.2, data=Location)) #79.45476 gam with a smooth term to year, two level factors with a plus after year
plot.gam(gam(TP.GLU~s(Year)+Location.2, data=Location))
fit.TP <- gam(TP.GLU~s(Year)+Location.2, data=Location)
AIC(gam(TP.GLU~s(Year), data=Location)) # 82.81808 gam with a smooth term to year, two level factors with a plus after year
AIC(gam(TP.GLU~s(Year, by=Location.2), data=Location)) #82.85721 gam with a smooth term to year, two level factors with a plus after year

########################     Time sereis to check seasonality and size       ############################

pdf(file="Results/Figures/MonthAnalysis.TPMean.pdf", width=6, height=5)
fit.1 <- gam(TPmean~s(Month, bs="cc", k=12), data=data)
b <- getViz(fit.1)
phe.m <- plot(b) +  l_points(shape = 19, size = 2.5, alpha = 0.2)+ l_fitLine(linetype = 3) +
  ggtitle(expression(paste(delta^15, "Trophic Position")))+
  l_ciLine(colour = 'red') + theme_classic() + labs(title = NULL)
print(phe.m, pages = 1)

dev.off()



pdf(file="Results/Figures/TPbyLength.pdf", width=8, height=4)

par(mfrow=c(1,1), mar=c(5,5,2,2))
length <- subset(data, Length>=1)
length(length$Sample.ID) #122
length.stand <- ifelse(length$Length>900, length$Length/10, length$Length)
length.st<- cbind(length, length.stand =length.stand)
plot(length.st$length.stand, length.st$TPmean,ylab="Trophic Position", xlab="Length (cm)", pch=16,
     ylim=c(2,6))
lm(length.st$TPmean~length.st$length.stand)
summary(lm(length.st$TPmean~length.st$length.stand))
abline(4.19, 0.000426, col='red')
dev.off()



########################     Creating Hierarchical Dataset       ############################
transform_to_log_scale <- function(x){
  if(x!=0){
    y <- (sign(x)) * (log(abs(x)))
  } else {
    y <-0
  }
  y 
}

transform_to_log <- function(x){
  
    y <- (sign(x)) * (log(abs(x)))
  y 
}

data.hGLU = data[,c("years","Location.2","TP.GLU", "Sample.ID")]
#data.hGLU =cbind(data.hGLU, AA=rep('Glu',length(data.hGLU$TP.GLU)), TP.norm=data.hGLU$TP.GLU-mean(na.omit(data.hGLU$TP.GLU)))
data.hGLU =cbind(data.hGLU, AA=rep('Glu',length(data.hGLU$TP.GLU)), TP.norm=((data.hGLU$TP.GLU)-mean(na.omit(data.hGLU$TP.GLU)))/sd(na.omit(data.hGLU$TP.GLU)))
data.hGLU =dplyr::rename(data.hGLU, TP = TP.GLU)

data.hASP = data[,c("years","Location.2","TP.ASP", "Sample.ID")]
#data.hASP =cbind(data.hASP, AA=rep('Asp',length(data.hASP$TP.ASP)), TP.norm=data.hASP$TP.ASP-mean(na.omit(data.hASP$TP.ASP)))
data.hASP =cbind(data.hASP, AA=rep('ASP',length(data.hASP$TP.ASP)), TP.norm=((data.hASP$TP.ASP)-mean(na.omit(data.hASP$TP.ASP)))/sd(na.omit(data.hASP$TP.ASP)))
data.hASP =dplyr::rename(data.hASP, TP = TP.ASP)

data.hPRO = data[,c("years","Location.2","TP.PRO", "Sample.ID")]
#data.hPRO =cbind(data.hPRO, AA=rep('Pro',length(data.hPRO$TP.PRO)), TP.norm=data.hPRO$TP.PRO-mean(na.omit(data.hPRO$TP.PRO)))
data.hPRO =cbind(data.hPRO, AA=rep('PRO',length(data.hPRO$TP.PRO)), TP.norm=((data.hPRO$TP.PRO)-mean(na.omit(data.hPRO$TP.PRO)))/sd(na.omit(data.hPRO$TP.PRO)))
data.hPRO =dplyr::rename(data.hPRO, TP = TP.PRO)

data.hVAL = data[,c("years","Location.2","TP.VAL", "Sample.ID")]
#data.hVAL =cbind(data.hVAL, AA=rep('Val',length(data.hVAL$TP.VAL)), TP.norm=data.hVAL$TP.VAL-mean(na.omit(data.hVAL$TP.VAL)))
data.hVAL =cbind(data.hVAL, AA=rep('VAL',length(data.hVAL$TP.VAL)), TP.norm=((data.hVAL$TP.VAL)-mean(na.omit(data.hVAL$TP.VAL)))/sd(na.omit(data.hVAL$TP.VAL)))
data.hVAL =dplyr::rename(data.hVAL, TP = TP.VAL)

data.hALA = data[,c("years","Location.2","TP.ALA", "Sample.ID")]
#data.hALA =cbind(data.hALA, AA=rep('Ala',length(data.hALA$TP.ALA)), TP.norm=data.hALA$TP.ALA-mean(na.omit(data.hALA$TP.ALA)))
data.hALA =cbind(data.hALA, AA=rep('ALA',length(data.hALA$TP.ALA)), TP.norm=((data.hALA$TP.ALA)-mean(na.omit(data.hALA$TP.ALA)))/sd(na.omit(data.hALA$TP.ALA)))
data.hALA =dplyr::rename(data.hALA, TP = TP.ALA)

data.hPHE = data[,c("years","PHE.mean", "Sample.ID")]
#data.hPHE =cbind(data.hPHE,PHE.norm=data.hPHE$PHE.mean-mean(na.omit(data.hPHE$PHE.mean)))
data.hPHE =cbind(data.hPHE,PHE.norm=(data.hPHE$PHE.mean-mean(na.omit(data.hPHE$PHE.mean)))/sd(na.omit(data.hPHE$PHE.mean)))


data.d13C = data[,c("years","d13C.s", "Sample.ID")]
#data.d13C =cbind(data.d13C,d13C.norm=data.d13C$d13C.s-mean(na.omit(data.d13C$d13C.s)))
data.d13C =cbind(data.d13C,d13C.norm=(data.d13C$d13C.s-mean(na.omit(data.d13C$d13C.s)))/sd(na.omit(data.d13C$d13C.s)))

########################     basic GAMS      ############################

#lets relationship be non linear, 
data.GLU.WA<- subset(data.hGLU, Location.2=='Inland'|Location.2=="Coastal")
summary(gam(TP~s(years, by=Location.2, k=3), data=data.GLU.WA)) #gam with a smooth term to year, two level factors with a plus after year
fit.1 <- gam(TP~Location.2+s(years, by=Location.2, k=3), data=data.GLU.WA) #gam with a smooth term to year, two level factors with a plus after year
AIC(fit.1)
AIC(gam(TP~s(years, by=Location.2, k=4), data=data.GLU.WA))

b<-getViz(fit.1)
pdf(file="Results/Presentation Figures/SmCoastDRAFT.pdf", width=5, height=5)
plot(sm(b,1)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ='#CCA65A',lwd=1) + 
  l_ciLine(colour = '#CCA65A', linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
dev.off()

pdf(file="Results/Presentation Figures/SmSalishDRAFT.pdf", width=5, height=5)
plot(sm(b,2)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ='#CCA65A',lwd=1) + 
  l_ciLine(colour = '#CCA65A', linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
dev.off()

pdf(file="Results/Presentation Figures/IntDRAFT.pdf", width=5, height=5)
plot(pterm(b, 1)) + l_ciBar(colour = c('#6FB1E7','#CCA65A'), lwd=1) + 
  l_fitPoints(colour=c("skyblue4",'#CCA65A'), cex=3)+
  l_points(shape = 19, size = 1, alpha = 0.1)
dev.off()

########################     Hierarchical Climate Models       ############################

data.hier <- rbind(data.hGLU, data.hASP)
data.hier <- rbind(data.hier, data.hPRO)
data.hier <- rbind(data.hier, data.hVAL)
data.hier <- rbind(data.hier, data.hALA)
data.hier<- subset(data.hier, Location.2=="Coastal"|Location.2=="Inland")
lag <- 1
data.hier <- cbind(data.hier, Year=data.hier$years+lag) # adding 1 year lag

dfa <- read.csv("Data/Compiled/DFA.csv")
dfa = dfa[,c("Year","Up.DFA1","Climate.DFA1.x", "SST.DFA1", "Dis.DFA1")]
data2<-merge(dfa,data.hier, by='Year', all.x=TRUE)
data2<-merge(data2,data.hPHE, by='Sample.ID', all.x=TRUE)
data2<-merge(data2,data.d13C, by='Sample.ID', all.x=TRUE)
write.csv(data2, 'Data/Compiled/Env.norm.csv')


ModelSelection.WA <- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~Climate.DFA1.x+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~SST.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Up.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Dis.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Climate.DFA1.x+SST.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Climate.DFA1.x+Dis.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~Dis.DFA1+SST.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~Up.DFA1+SST.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~Up.DFA1+Dis.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AIC(lmer(y~Up.DFA1+SST.DFA1+Dis.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#7
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+Dis.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#8
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+SST.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#9
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+Dis.DFA1+SST.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#10
                      AIC(lmer(y~Climate.DFA1.x+SST.DFA1+Dis.DFA1+Location.2+(1|AA)+(1|Sample.ID), data=dataframe))
  )
  
 
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Climate", "SST", "Up", "Dis","1. Clim, SST", "2. Clim, Dis", "3. Clim, Up", "4. SST, Dis",
                   "5. SST, Up", "6. Up, Dis", "7. SST, Dis, Up", "8. Clim, Up, Dis", "9. Clim, Up, SST", "10. Clim, Up,
                   SST, Dis", "11. Clim Dis SST"
                )
  
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
model.selectionENV <- ModelSelection.WA(data2, n, data2$TP)
model.selectionENV
length(model.selectionENV[,1])



modelENV<- lmer(TP~Location.2+(1|AA), data=data2)
summary(modelENV)
summary(lmer(TP~Climate.DFA1.x+Location.2+(1|AA), data=data2))

coef(modelENV)
plot(modelENV)

subset()
########################     Hierarchical PREY Models       ############################

prey <- read.csv("Data/Compiled/WA.Prey2.tot.csv")
data.merge <- list(prey,data.hier, dfa)
data3 <- Reduce(function(x,y) merge(x,y, by='Year', all.x=TRUE), data.merge)
data3<-merge(data3,data.hPHE, by='Sample.ID', all.x=TRUE)
data3<-merge(data3,data.d13C, by='Sample.ID', all.x=TRUE)
data3 <- subset(data3, Year>=1973&Year<=2008)


ModelSelection.WA2 <- function(dataframe,n, y) {
  
  aic.output <- rbind(AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Herring.Biomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AIC(lmer(y~allSmolt+Chinook+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#7
                      AIC(lmer(y~Herring.Biomass+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#8
                      AIC(lmer(y~Herring.Biomass+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#9
                      AIC(lmer(y~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      
                      AIC(lmer(y~Location.2+Chum+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~Chum+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chum+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chum+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chum+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Chum+Herring.Biomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Chum+Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~Chum+Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~Chum+HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~Chum+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~Chum+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+Coho+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~Coho+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Coho+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Coho+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Coho+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~Coho+Herring.Biomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~Coho+Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~Coho+Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~Coho+HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~Coho+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~Coho+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~Location.2+HarborSeal+(1|AA)+(1|Sample.ID), data=dataframe)), 
                      AIC(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HarborSeal+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HarborSeal+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HarborSeal+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),
                      AIC(lmer(y~HarborSeal+Herring.Biomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
                      AIC(lmer(y~HarborSeal+Herring.Biomass+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
                      AIC(lmer(y~HarborSeal+Herring.Biomass+allSmolt+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
                      AIC(lmer(y~HarborSeal+HakeBiomass+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~HarborSeal+allSmolt+Chinook+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~HarborSeal+allSmolt+HakeBiomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      
                      AIC(lmer(y~HarborSeal+Chum+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
                      AIC(lmer(y~HarborSeal+Coho+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
                      AIC(lmer(y~Chum+Coho+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
                      AIC(lmer(y~Chum+Coho+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe))#6
                      
                      

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

model.selectionPREY <- ModelSelection.WA2(data3, n, data3$TP)
length(model.selectionPREY[,1])
model.selectionPREY
modelPREY<-lmer(TP~HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=data3)
summary(modelPREY)
coef(modelPREY)



########################     Hierarchical Nutrient Models       ############################


ModelSelection.WA3 <- function(dataframe,n, y) {
  
  aic.output <- rbind(
    AIC(lmer(TP~Location.2+(1|AA), data=dataframe)),
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
model.selectionNUTRIENT <- ModelSelection.WA3(data2, n, data2$TP)
model.selectionNUTRIENT
length(model.selectionNUTRIENT[,1])

modelNUTRIENT<-lmer(TP~PHE.norm+d13C.norm+Location.2+(1|AA)+(1|Sample.ID),data=data3)
summary(modelNUTRIENT)


########################     Hierarchical FULL Models       ############################


ModelSelection.WAFULL <- function(dataframe,n, y) {
  
  aic.output <- rbind(
    AIC(lmer(y~Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#1
    AIC(lmer(y~PHE.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#2
    AIC(lmer(y~PHE.norm+Climate.DFA1.x+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#3
    AIC(lmer(y~PHE.norm+d13C.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#4
    AIC(lmer(y~PHE.norm+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#5
    AIC(lmer(y~PHE.norm+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#6
    
    AIC(lmer(y~d13C.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#7
    AIC(lmer(y~d13C.norm+Climate.DFA1.x+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#8
    AIC(lmer(y~d13C.norm+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#9
    AIC(lmer(y~d13C.norm+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#10
    
    AIC(lmer(y~Climate.DFA1.x+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#11
    AIC(lmer(y~Climate.DFA1.x+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#12
    AIC(lmer(y~Climate.DFA1.x+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#13
    
    AIC(lmer(y~HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#14
    AIC(lmer(y~HarborSeal+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#15
    AIC(lmer(y~Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#16
    
    AIC(lmer(y~PHE.norm+Climate.DFA1.x+d13C.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#17
    AIC(lmer(y~PHE.norm+Climate.DFA1.x+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#18
    AIC(lmer(y~PHE.norm+Climate.DFA1.x+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#19
    AIC(lmer(y~PHE.norm+d13C.norm+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#20
    AIC(lmer(y~PHE.norm+Herring.Biomass+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#21
    
    AIC(lmer(y~Climate.DFA1.x+d13C.norm+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#22
    AIC(lmer(y~Climate.DFA1.x+d13C.norm+Herring.Biomass+Location.2+(1|AA)+(1|Sample.ID)+(1|Sample.ID), data=dataframe)),#23
    AIC(lmer(y~Climate.DFA1.x+Herring.Biomass+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#24
   
    AIC(lmer(y~Herring.Biomass+HarborSeal+d13C.s+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#25
    
    AIC(lmer(y~Herring.Biomass+HarborSeal+d13C.s+PHE.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#26
    AIC(lmer(y~Herring.Biomass+HarborSeal+d13C.s+Climate.DFA1.x+Location.2+(1|AA)+(1|Sample.ID), data=dataframe)),#27
    AIC(lmer(y~Herring.Biomass+HarborSeal+d13C.s+Climate.DFA1.x+PHE.norm+Location.2+(1|AA)+(1|Sample.ID), data=dataframe))#28
    
  )
  
  names <- seq(1,n,1)
  model.names <- c( "1. Location","2. Phe", "3. PHE Climate", "4. PHE, d13C", "5. PHE, HarborSeal",
                    "6. PHE, Herring", "7. d13C", "8. d13C, Climate", "9. d13C, HarborSeal",
                    "10. d13C, Herring", "11. Climate", "12. Climate, Harbor Seal", "13. Climate, Herring",
                    "14.Harbor Seal", "15. Harbor Seal, Herring", "16. Herring", "17. PHE, Climate, d13C",
                    "18. PHE, Climate, Herring", "19. PHE, Climate, Herring", "20. PHE, d13C, Harbor Seal",
                    "21. PHE, Herring, HarborSeal", "22. Climate, d13C, Harbor Seal", "23. Climate, d13C, Herring",
                    "24. Climate, Herring, Harbor Seal", "25. Herring, Harbor Seal, d13C", "26. Herring, Seal, d13C, PHE",
                    "27. Herring, Seal, d13C, Climate", "28. Herring, Seal, Climate, PHE")
  
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
model.selectionFULL <- ModelSelection.WAFULL(data3, n, data3$TP)
model.selectionFULL



modelFULL<-lmer(TP~PHE.norm+d13C.norm+HarborSeal+Location.2+(1|AA)+(1|Sample.ID), data=data3)
vif(modelFULL)

summary(modelFULL)
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1|AA),data=data3))
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1+Location.2|AA),data=data3))
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1+Location.2+PHE.norm|AA),data=data3))
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1+Location.2+d13C.norm|AA),data=data3))
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1+Location.2+HarborSeal|AA),data=data3))
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1+Location.2+PHE.norm+HarborSeal|AA),data=data3))
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1+Location.2+PHE.norm+d13C.norm|AA),data=data3))
summary(lmer(TP~PHE.norm+d13C.norm+Location.2+(1+Location.2+PHE.norm+d13C.norm+HarborSeal|AA),data=data3))



coef(modelNUTRIENT)
#AIC(lmer(TP.norm~PHE.mean+(1|AA), data=dataframe)),#10

########################    VIF
vif(modelNUTRIENT)
vif(modelPREY)
vif(modelENV)
vif(modelFULL)




########################    Mechanistic Models       ############################


ModelSelection.WAHake <- function(dataframe,n) {
  
  aic.output <- rbind(
    AIC(lm(HakeBiomass~PHE.norm, data=dataframe)),
    AIC(lm(HakeBiomass~Climate.DFA1.x, data=dataframe)),
    AIC(lm(HakeBiomass~Climate.DFA1.x+PHE.norm, data=dataframe)),
    AIC(lm(HakeBiomass~Climate.DFA1.x*PHE.norm, data=dataframe)),
    
    AIC(lm(HakeBiomass~TP.norm, data=dataframe)),
    AIC(lm(HakeBiomass~TP.norm+PHE.norm, data=dataframe)),
    AIC(lm(HakeBiomass~TP.norm+Climate.DFA1.x, data=dataframe)),
    
    AIC(lm(HakeBiomass~HarborSeal, data=dataframe)),
    AIC(lm(HakeBiomass~HarborSeal+PHE.norm, data=dataframe)),
    AIC(lm(HakeBiomass~HarborSeal+Climate.DFA1.x, data=dataframe))
  )
  
  names <- seq(1,n,1)
  model.names <- c( "PHE", "Climate", "Climate, PHE","Climate*PHE", "TP", "TP, PHE", "TP, Climate", "Seal", "Seal, PHE", "Seal, Climate")
  
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
model.selectionHAKE <- ModelSelection.WAHake(data4, n)
model.selectionHAKE

ModelSelection.WAHerring <- function(dataframe,n) {
  
  aic.output <- rbind(
    AIC(lm(Herring.Biomass~PHE.norm, data=dataframe)),
    AIC(lm(Herring.Biomass~Climate.DFA1.x, data=dataframe)),
    AIC(lm(Herring.Biomass~Climate.DFA1.x+PHE.norm, data=dataframe)),
    AIC(lm(Herring.Biomass~Climate.DFA1.x*PHE.norm, data=dataframe)),
    
    AIC(lm(Herring.Biomass~d13C.norm, data=dataframe)),
    AIC(lm(Herring.Biomass~Climate.DFA1.x+d13C.norm, data=dataframe)),
    AIC(lm(Herring.Biomass~d13C.norm+Climate.DFA1.x*PHE.norm, data=dataframe)),
    
    AIC(lm(Herring.Biomass~TP.norm, data=dataframe)),
    AIC(lm(Herring.Biomass~TP.norm+PHE.norm, data=dataframe)),
    AIC(lm(Herring.Biomass~TP.norm+Climate.DFA1.x, data=dataframe)),
    AIC(lm(Herring.Biomass~TP.norm+PHE.norm+Climate.DFA1.x, data=dataframe)),
    AIC(lm(Herring.Biomass~HarborSeal, data=dataframe)),
    AIC(lm(Herring.Biomass~HarborSeal+PHE.norm, data=dataframe)),
    AIC(lm(Herring.Biomass~HarborSeal+Climate.DFA1.x+HakeBiomass, data=dataframe)),
    AIC(lm(Herring.Biomass~HarborSeal+Climate.DFA1.x*PHE.norm+HakeBiomass, data=dataframe)),
    AIC(lm(Herring.Biomass~HarborSeal+Climate.DFA1.x+PHE.norm+HakeBiomass, data=dataframe))
    
  )
  
  names <- seq(1,n,1)
  model.names <- c( "PHE", "Climate", "Climate, PHE", "PHE*Climate", "d13C", "Climate, d13C", 
                    "climate*d13c", "TP", "TP, PHE", "TP, Climate", "TP, Climate, PHE", "Seal", "Seal, PHE", "Seal, Climate",
                    "Seal, climate*phe", "Seal, climate, phe")
  
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
model.selectionHERRING <- ModelSelection.WAHerring(data4, n)
model.selectionHERRING

summary(lm(HakeBiomass~PHE.norm, data=data4))
summary(lm(HakeBiomass~Climate.DFA1.x, data=data4))
summary(lm(HakeBiomass~Climate.DFA1.x*PHE.norm, data=data4))


summary(lm(Herring.Biomass~PHE.norm, data=data4))
summary(lm(Herring.Biomass~Climate.DFA1.x, data=data4))
summary(lm(Herring.Biomass~Climate.DFA1.x*PHE.norm, data=data4))
summary(lm(Herring.Biomass~HakeBiomass, data=data4))



ModelSelection.WA4 <- function(dataframe,n) {
  
  aic.output <- rbind(AIC(lmer(TP.norm~Location.2+(1|AA), data=dataframe)), #1
                      AIC(lmer(TP.norm~PHE.norm+Location.2+(1|AA), data=dataframe)), #2
                      AIC(lmer(TP.norm~PHE.norm+Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(TP.norm~PHE.norm+HakeBiomass+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(TP.norm~PHE.norm+Herring.Biomass+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(TP.norm~PHE.norm+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(TP.norm~PHE.norm+Herring.Biomass+Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),#7
                      AIC(lmer(TP.norm~PHE.norm+Climate.DFA1.x+HakeBiomass+Location.2+(1|AA), data=dataframe)),#8
                      
                      AIC(lmer(TP.norm~Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),#9
                      AIC(lmer(TP.norm~Climate.DFA1.x+HakeBiomass+Location.2+(1|AA), data=dataframe)),#10
                      AIC(lmer(TP.norm~Climate.DFA1.x+Herring.Biomass+Location.2+(1|AA), data=dataframe)),#11
                      AIC(lmer(TP.norm~Climate.DFA1.x+Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#12
                      
                      AIC(lmer(TP.norm~HakeBiomass+Location.2+(1|AA), data=dataframe)),#13
                      AIC(lmer(TP.norm~Herring.Biomass+Location.2+(1|AA), data=dataframe)),#14
                      AIC(lmer(TP.norm~Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#15
                      AIC(lmer(TP.norm~Climate.DFA1.x+Herring.Biomass+HakeBiomass+PHE.norm+Location.2+(1|AA), data=dataframe))#16
                      
                      
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("1. Location", "2. Phe", "3. Phe, Climate", "4. Phe, Hake", "5. Phe, Herring", "6. Phe, Hake, Herring", 
                   "7. PHE, Climate, Herring", "8. Phe, Climate, Hake", "9. Climate", "10. Climate, Hake", "11. Climate, Herring",
                   "12. Climate, Herring, Hake", "13. Hake", "14. Herring", "15. Herring, Hake", "16. Herring, Climate, Phe, Hake")
  
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
model.selectionFULL <- ModelSelection.WA4(data3, n)
model.selectionFULL

########################     Lasso Regression       ############################
library(dplyr)
library(lme4)

data5<-subset(data3, AA=="Glu")

data4 <- data5 %>%
  select(-TP) %>%
  select(-Sample.ID) %>%
  select(-Year) %>%
  select(-X) %>%
  select(-years.x) %>%
  select(-years.y) %>%
  select(-AA)%>%
  select(-EulachonLandings)%>%
  select(-WildProduction)%>%
  select(-HatcherySmolts)

Wa.lasso <- data4[complete.cases(data4), ]


y<- Wa.lasso$TP.norm
Wa.lasso <- Wa.lasso %>%
  select(-PHE.mean)%>%
  select(-TP.norm)%>%
  select(-d13C.s)%>%
  select(-years)%>%
  #select(-PHE.norm)%>%
  #select(-HarborSeal)%>%
 # select(-Climate.DFA1.x)%>%
 select(-Up.DFA1)%>%
 select(-SST.DFA1)%>%
  select(-Dis.DFA1)
 # select(-d13C.norm)

i_num <- sapply(Wa.lasso, is.numeric)
Wa.lasso[, i_num] <- apply(Wa.lasso[, i_num], 2, function(x) ifelse(is.na(x), median(x, na.rm=T), x))

# impute empty categories
Wa.lasso[, !i_num] <- apply(Wa.lasso[, !i_num, drop=F], 2, function(x) {
  x[x==""] <- "empty"
  x[is.na(x)] <- "missing"
  x
})
#nlevels(X$Location)

# get the numLevels vector containing the number of categories
X <- Wa.lasso


numLevels <- X %>% sapply(nlevels)
numLevels[8] <- 2
#apply(X[, 'Location'], factor)
#numLevels[9] <- 2
numLevels[numLevels==0] <- 1


#X[8][X[8] == 2] <- 0
#X[8][X[8] == 3] <- 1

# make the categorical variables take integer values starting from 0
X[, !i_num] <- apply(X[, !i_num], 2, function(col) as.integer(as.factor(col)) - 1)
sapply(X, class)






library(glinternet)
set.seed(1001)

head(X)

pairs <- matrix(c(8,1, 8,2, 8,3, 8,4,
                  8,5, 8,6, 8,7), ncol=7, byrow=TRUE)

cv_fit <- glinternet.cv(X, y, numLevels=numLevels, interactionPairs = pairs)
cv_fit <- glinternet.cv(X, y, numLevels=numLevels)

plot(cv_fit)
i_1Std <- which(cv_fit$lambdaHat1Std == cv_fit$lambda)
coefs <- coef(cv_fit$glinternetFit)[[i_1Std]]


coefs$mainEffects
idx_num <- (1:length(i_num))[i_num]
idx_cat <- (1:length(i_num))[!i_num]
names(numLevels)[idx_cat[coefs$mainEffects$cat]]
names(numLevels)[idx_num[coefs$mainEffects$cont]]
coefs$mainEffects
coefs$mainEffectsCoef

coefs$interactions
coefs$interactionsCoef




-0.8836366+10.89 #coastal intercept
10.89+0.4197787 #salish sea intercept





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

