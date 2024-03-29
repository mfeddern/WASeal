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
data2 <- read.csv("Data/Compiled/WASealAAandTP2.csv")
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

data2 <- read.csv("Data/Compiled/WASealAAandTP2.csv")
data <- read.csv("Data/Compiled/SI Data.csv")


########################     Calculating TEFs                      ############################

ger <- read.csv("Data/Compiled/GermainData2.csv")
ger.sealAA <- ger[1:7,7:19]
#TEF.JN<- read.csv("Data/Compiled/TEFCh.csv")
JN<- read.csv("Data/Compiled/Nielsen_Beta.csv")
Beta.JN <- JN[1:2,]
TEF.JN <- JN[3:4,]


find.numeric <- sapply(ger.sealAA, is.numeric)
meanAA <- colMeans(ger.sealAA[, find.numeric])
herring <- ger[8,7:19]
find.numeric <- sapply(herring, is.numeric)
herring <- colMeans(herring[, find.numeric])
phe<- herring['Phe']
TEF <- (meanAA-rep(meanAA[10],length(meanAA)))-herring+rep(herring[10],length(herring))
herringtp <- (((herring['Glu']-herring['Phe'])-Beta.JN$GLU[1]-TEF['Glu'])/TEF.JN$GLU[1])+2
herringtp.VAL <- (((herring['Val']-herring['Phe'])-Beta.JN$VAL [1]-TEF['Val'])/TEF.JN$VAL[1])+2
herringtp.Pro <- (((herring['Pro']-herring['Phe'])-Beta.JN$PRO [1]-TEF['Pro'])/TEF.JN$PRO[1])+2
herringtp.Asp <- (((herring['Asp']-herring['Phe'])-Beta.JN$ASP [1]-TEF['Asp'])/TEF.JN$ASP[1])+2
herringtp.Ala <- (((herring['Ala']-herring['Phe'])-Beta.JN$ALA [1]-TEF['Ala'])/TEF.JN$ALA[1])+2

HerringTP<-cbind(herringtp, herringtp.VAL, herringtp.Pro, herringtp.Asp, herringtp.Ala)
colnames(HerringTP)<- c("Glu", "Val", "Pro","Asp", "Ala")





#write.csv(data4, 'Data/Compiled/HierarchicalData.csv')

TP.GLU <- ((data$GLU.mean-data$PHE.mean-Beta.JN$GLU[1]-TEF['Glu'])/TEF.JN$GLU[1])+2
plot(TP.GLU, ylim = c(0,6))
mean(na.omit(TP.GLU))


TP.ALA <- (((data$ALA.mean-data$PHE.mean)-Beta.JN$ALA[1]-TEF['Ala'])/TEF.JN$ALA[1])+2
mean(na.omit(TP.ALA))

TP.PRO <- (((data$PRO.mean-data$PHE.mean)-Beta.JN$PRO[1]-TEF['Pro'])/TEF.JN$PRO[1])+2
mean(na.omit(TP.PRO))

TP.ASP <- (((data$ASP.mean-data$PHE.mean)-Beta.JN$ASP[1]-TEF.JN$ASP[1])/TEF['Asp'])+2
plot(TP.ASP, ylim = c(0,6))


TP.VAL <- (((data$VAL.mean-data$PHE.mean)-Beta.JN$VAL[1]-TEF['Val'])/TEF.JN$VAL[1])+2
mean(na.omit(TP.VAL))

data<-cbind(data, TP.ALA, TP.ASP, TP.GLU, TP.PRO,
            TP.VAL)

########################     Summary Statistics                      ############################


par(mfrow=c(1,1), oma=c(0,0,0,0))
plot(data$years, data$TP.GLU, pch=16, ylab="Trophic Position", xlab="Year",ylim=c(2,5.5))
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



plot(In$years, In$TP.GLU, pch=16, col='blue',ylab="Trophic Position",xlab="Year", ylim=c(2,5.5))
points(Co$years, Co$TP.GLU, pch=16, col='green3')
legend(1925, 5.63, legend=c("Salish Sea", "Coastal"),
       col=c("blue", "green3"), pch=16, cex=0.8)
summary(lm(TP.GLU~years+Location.2, data=Location))
t.test(TP.GLU~Location.2, data=Location) #p=0.01252




summary(gam(TP.GLU~s(Year)+Location.2, data=Location)) #gam with a smooth term to year, two level factors with a plus after year
AIC(gam(TP.GLU~s(Year)+Location.2, data=Location)) #79.45476 gam with a smooth term to year, two level factors with a plus after year
plot.gam(gam(TP.GLU~s(Year)+Location.2, data=Location))
fit.TP <- gam(TP.GLU~s(Year)+Location.2, data=Location)
AIC(gam(TP.GLU~s(Year), data=Location)) # 82.81808 gam with a smooth term to year, two level factors with a plus after year
AIC(gam(TP.GLU~s(Year, by=Location.2), data=Location)) #82.85721 gam with a smooth term to year, two level factors with a plus after year

########################     Time sereis to check seasonality and size       ############################

pdf(file="Results/Figures/MonthAnalysis.TPMean.pdf", width=6, height=5)
fit.1 <- gam(TP.ger~s(Month, bs="cc", k=12), data=data2)
b <- getViz(fit.1)
phe.m <- plot(b) +  l_points(shape = 19, size = 2.5, alpha = 0.2)+ l_fitLine(linetype = 3) +
  ggtitle(expression(paste(delta^15, "Trophic Position")))+
  l_ciLine(colour = 'red') + theme_classic() + labs(title = NULL)
print(phe.m, pages = 1)

dev.off()

palette(c('black','#6FB1E7','#D494E1','#CCA65A','#7EBA68','#00C1B2'))


pdf(file="Results/Figures/TPbyLength.pdf", width=8, height=4)
shapes = c(16, 17, 18) 
shapes <- shapes[as.numeric(length$Location.2)]
par(mfrow=c(1,1), mar=c(5,5,2,2))
length <- subset(data, Length>=1&Location.2=="Inland"|Location.2=="Coastal")
length(length$Sample.ID) #73
length.stand <- ifelse(length$Length>900, length$Length/10, length$Length)
length.st<- cbind(length, length.stand =length.stand)
plot(log(length.st$length.stand), length.st$TP.GLU,ylab="Trophic Position", xlab="Length (cm)", 
     pch=shapes,
     ylim=c(2,6), col=length$Sex, bty='n')
lm(length.st$TP.GLU~log(length.st$length.stand))
summary(lm(length.st$TP.GLU~log(length.st$length.stand)))
abline(3.91, 0.00149, col='red')
dev.off()

fit.2 <- gam(TP.GLU~s(length.stand, k=3), data=length.st) #82.85721 gam with a smooth term to year, two level factors with a plus after year
b <- getViz(fit.2)
length <- plot(b) +  l_points(shape = 19, size = 2.5, alpha = 0.2)+ l_fitLine(linetype = 3) +
  ggtitle(expression(paste(delta^15, "Trophic Position")))+
  l_ciLine(colour = 'red') + theme_classic() + labs(title = NULL)
summary(fit.2)
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

data.hGLU = data[,c("years","Location.2","TP.GLU", "Sample.ID", "Sex", "Length")]
#data.hGLU =cbind(data.hGLU, AA=rep('Glu',length(data.hGLU$TP.GLU)), TP.norm=data.hGLU$TP.GLU-mean(na.omit(data.hGLU$TP.GLU)))
data.hGLU =cbind(data.hGLU, AA=rep('Glu',length(data.hGLU$TP.GLU)), TP.norm=((data.hGLU$TP.GLU)-mean(na.omit(data.hGLU$TP.GLU)))/sd(na.omit(data.hGLU$TP.GLU)))
data.hGLU =dplyr::rename(data.hGLU, TP = TP.GLU)

data.hASP = data[,c("years","Location.2","TP.ASP", "Sample.ID","Sex", "Length")]
#data.hASP =cbind(data.hASP, AA=rep('Asp',length(data.hASP$TP.ASP)), TP.norm=data.hASP$TP.ASP-mean(na.omit(data.hASP$TP.ASP)))
data.hASP =cbind(data.hASP, AA=rep('ASP',length(data.hASP$TP.ASP)), TP.norm=((data.hASP$TP.ASP)-mean(na.omit(data.hASP$TP.ASP)))/sd(na.omit(data.hASP$TP.ASP)))
data.hASP =dplyr::rename(data.hASP, TP = TP.ASP)

data.hPRO = data[,c("years","Location.2","TP.PRO", "Sample.ID","Sex", "Length")]
#data.hPRO =cbind(data.hPRO, AA=rep('Pro',length(data.hPRO$TP.PRO)), TP.norm=data.hPRO$TP.PRO-mean(na.omit(data.hPRO$TP.PRO)))
data.hPRO =cbind(data.hPRO, AA=rep('PRO',length(data.hPRO$TP.PRO)), TP.norm=((data.hPRO$TP.PRO)-mean(na.omit(data.hPRO$TP.PRO)))/sd(na.omit(data.hPRO$TP.PRO)))
data.hPRO =dplyr::rename(data.hPRO, TP = TP.PRO)

data.hVAL = data[,c("years","Location.2","TP.VAL", "Sample.ID", "Sex", "Length")]
#data.hVAL =cbind(data.hVAL, AA=rep('Val',length(data.hVAL$TP.VAL)), TP.norm=data.hVAL$TP.VAL-mean(na.omit(data.hVAL$TP.VAL)))
data.hVAL =cbind(data.hVAL, AA=rep('VAL',length(data.hVAL$TP.VAL)), TP.norm=((data.hVAL$TP.VAL)-mean(na.omit(data.hVAL$TP.VAL)))/sd(na.omit(data.hVAL$TP.VAL)))
data.hVAL =dplyr::rename(data.hVAL, TP = TP.VAL)

data.hALA = data[,c("years","Location.2","TP.ALA", "Sample.ID","Sex", "Length")]
#data.hALA =cbind(data.hALA, AA=rep('Ala',length(data.hALA$TP.ALA)), TP.norm=data.hALA$TP.ALA-mean(na.omit(data.hALA$TP.ALA)))
data.hALA =cbind(data.hALA, AA=rep('ALA',length(data.hALA$TP.ALA)), TP.norm=((data.hALA$TP.ALA)-mean(na.omit(data.hALA$TP.ALA)))/sd(na.omit(data.hALA$TP.ALA)))
data.hALA =dplyr::rename(data.hALA, TP = TP.ALA)

data.hPHE = data[,c("years","PHE.mean", "Sample.ID","Sex", "Length")]
#data.hPHE =cbind(data.hPHE,PHE.norm=data.hPHE$PHE.mean-mean(na.omit(data.hPHE$PHE.mean)))
data.hPHE =cbind(data.hPHE,PHE.norm=(data.hPHE$PHE.mean-mean(na.omit(data.hPHE$PHE.mean)))/sd(na.omit(data.hPHE$PHE.mean)))


data.d13C = data[,c("years","d13C.s", "Sample.ID","Sex", "Length")]
#data.d13C =cbind(data.d13C,d13C.norm=data.d13C$d13C.s-mean(na.omit(data.d13C$d13C.s)))
data.d13C =cbind(data.d13C,d13C.norm=(data.d13C$d13C.s-mean(na.omit(data.d13C$d13C.s)))/sd(na.omit(data.d13C$d13C.s)))



data.hier <- rbind(data.hGLU, data.hASP)
data.hier <- rbind(data.hier, data.hPRO)
data.hier <- rbind(data.hier, data.hVAL)
data.hier <- rbind(data.hier, data.hALA)
data.hier<- subset(data.hier, Location.2=="Coastal"|Location.2=="Inland")
lag <- 1
data.hier <-cbind(Year=data.hier$years+lag, data.hier)
total3 <- left_join(data.hier, data.d13C, by="Sample.ID")
total3 <- left_join(total3,data.hPHE, by="Sample.ID")

prey <- read.csv("Data/Compiled/WA.Prey2.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(total3, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData_2.csv')

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



########################     Germain AAs      ############################



TP.GLU.g <- (((ger.sealAA$Glu-ger.sealAA$Phe)-Beta.JN$GLU[1])/4.3)+1
plot(TP.GLU.g, ylim = c(0,6))
mean(TP.GLU.g)

TP.ALA.g <- (((ger.sealAA$Glu-ger.sealAA$Phe)-Beta.JN$ALA[1])/3.33)+1
plot(TP.ALA.g, ylim = c(0,10))
mean(TP.ALA.g)

#TP.ILE <- (((data$GLU.mean-data$PHE.mean)-Beta.JN$ILE-TEF$Glu)/6.6)+2
#plot(TP.GLU, ylim = c(0,6))

TP.PRO.g <- (((ger.sealAA$Pro-ger.sealAA$Phe)-Beta.JN$PRO[1])/6.68)+1
plot(TP.PRO.g, ylim = c(0,6))
mean(TP.PRO.g)

#TP.LEU <- (((data$LEU.mean-data$PHE.mean)-Beta.JN$LEU[1])/TEF$Leu)+1
#plot(TP.LEU, ylim = c(0,6))

TP.ASP.g <- (((ger.sealAA$Asp-ger.sealAA$Phe)-Beta.JN$ASP[1])/4.08)+1
plot(TP.ASP.g, ylim = c(0,6))
mean(TP.ASP.g)

TP.VAL.g <- (((ger.sealAA$Val-ger.sealAA$Phe)-Beta.JN$VAL[1])/8.31)+1
plot(TP.VAL.g, ylim = c(0,6))
mean(TP.VAL.g)

data2<-cbind(data, TP.ALA, TP.ASP, TP.GLU, TP.PRO,
            TP.VAL)


############### Time Series #################


dataTS<-data2 %>% select(
  TP.ALA,
  TP.GLU,
  TP.PRO,
  TP.VAL,
  Location.2,
  Year)





dataTS.d13C <- dataTS.d13C[complete.cases(dataTS.d13C), ]
dataTS.d13C<- subset(dataTS.d13C, d15N>10)
dataTS.d13C$Location<- factor(dataTS.d13C$Location, levels = c("a.BB","b.SC",
                                                               "c.SE","d.SS","e.C"),
                              labels = c("Eastern Bering Sea", "Northern GoA", "Southeast GoA", "Salish Sea", "Coastal WA")
)
palette(c("#ED90A4", "#C0AB52","#28BBD7", "#4FBF85",  "#C699E7"))
palette()
d13C.TS <- ggplot(dataTS.d13C, aes(x = Year, y = d13C)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  geom_point(aes(color = Location, alpha=0.5), pch=16, size=2.5) +
  facet_grid(Location~.,labeller = labeller(dataTS.d13C$Location.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)


d15N.TS<-ggplot(dataTS.d13C, aes(x = Year, y = d15N)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  geom_point(aes(color = Location, alpha=0.5), pch=16, size=2.5) +
  facet_grid(Location~.,labeller = labeller(dataTS.d13C$Location.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)





dataTS.phe <-data %>% select(
  PHE.mean,
  Location.2,
  Year)%>%
  mutate(Location =recode(Location.2, 'BB' = "a.BB",'SC'="b.SC",
                          'SE'="c.SE", 'Inland'="d.SS", 'Coastal'="e.C"))

dataTS.phe <- dataTS.phe[complete.cases(dataTS.phe), ]

dataTS.phe$Location<- factor(dataTS.phe$Location, levels = c("a.BB","b.SC",
                                                             "c.SE","d.SS","e.C"),
                             labels = c("Eastern Bering Sea", "Northern GoA", "Southeast GoA", "Salish Sea", "Coastal WA")
)
palette(c("#ED90A4", "#C0AB52","#28BBD7", "#4FBF85",  "#C699E7"))
palette()
PHE.TS <- ggplot(dataTS.phe, aes(x = Year, y = PHE.mean)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  geom_point(aes(color = Location, alpha=0.5), pch=16, size=2.5) +
  facet_grid(Location~.,labeller = labeller(dataTS.d13C$Location.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)


PHE.TS <- ggplot(dataTS.phe, aes(x = Year, y = PHE.mean)) + 
  #stat_smooth(fullrange=TRUE,aes(color = Location.2)) +
  geom_point(aes(color = Location, alpha=0.5), pch=16, size=2.5) +
  # facet_grid(Location~.,labeller = labeller(dataTS.d13C$Location.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  guides(colour=FALSE, alpha=FALSE)


pdf(file="Results/Figures/TS.pheRAW.pdf", width=5, height=6.5)
PHE.TS
dev.off()
pdf(file="Results/Figures/TS.d13cRAW.ns.pdf", width=5, height=6.5)
d13C.TS
dev.off()
pdf(file="Results/Figures/TS.d15NRAW.pdf", width=5, height=6.5)
d15N.TS
dev.off()

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


