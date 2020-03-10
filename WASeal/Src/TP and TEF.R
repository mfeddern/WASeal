rm(list = ls())
library(mgcViz)
library(MARSS)
library(forecast)
library(EBImage)
library(dplyr)
library(lme4)
#install.packages("digitize")
#install.packages("mgcv")
library(mgcv)
data <- read.csv("Data/Compiled/WASealAAandTP2.csv")
#install.packages("visreg")
library(visreg)
#data <- read.csv("SortedTP.csv")
#dev.off()
#################################################################################################
########################     Calculating TEFs                      ############################
#################################################################################################

ger <- read.csv("Data/Compiled/GermainData2.csv")
ger.sealAA <- ger[1:7,7:19]

find.numeric <- sapply(ger.sealAA, is.numeric)
meanAA <- colMeans(ger.sealAA[, find.numeric])
herring <- ger[8,7:19]
find.numeric <- sapply(herring, is.numeric)
herring <- herring[, find.numeric]
PHE <- rep(meanAA[10]-herring[10], length(meanAA))
TEF <- (meanAA-herring) - PHE

TEF.plank <- read.csv("Data/Compiled/Mclelland.csv")
TP.GLU <- (((data$GLU.mean-data$PHE.mean)-4.3-3.4)/TEF$Glu)+2
TP.ALA <- (((data$GLU.mean-data$PHE.mean)-2.5-3.4)/TEF$Ala)+2
TP.ASP <- (((data$GLU.mean-data$PHE.mean)-3.5-3.4)/TEF$Asp)+2
TP.VAL <- (((data$GLU.mean-data$PHE.mean)-7.5-3.4)/TEF$Val)+2
TP.PRO <- (((data$GLU.mean-data$PHE.mean)-5.5-3.4)/TEF$Pro)+2

TP.AA <- cbind(TP.GLU, TP.ALA, TP.ASP, TP.VAL, TP.PRO)
TP.AA <- cbind(TP.AA, TPmean=rowMeans(TP.AA, na.rm = FALSE))
data <- cbind(data,TP.AA)

write.csv(data, 'TPData_nonhier.csv')
#################################################################################################
########################     Summary Statistics                      ############################
#################################################################################################


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

#lets relationship be non linear, 
summary(gam(TP.GLU~s(Year, by=Location.2), data=Location)) #gam with a smooth term to year, two level factors with a plus after year
#################################################################################################
########################     Time sereis to check seasonality and size       ############################
#################################################################################################

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


data.hGLU = data[,c("years","Location.2","TP.GLU")]
data.hGLU =cbind(data.hGLU, AA=rep('Glu',length(data.hGLU$TP.GLU)), TP.norm=transform_to_log_scale(data.hGLU$TP.GLU-mean(na.omit(data.hGLU$TP.GLU))))
data.hGLU =dplyr::rename(data.hGLU, TP = TP.GLU)

data.hASP = data[,c("years","Location.2","TP.ASP")]
data.hASP =cbind(data.hASP, AA=rep('Asp',length(data.hASP$TP.ASP)), TP.norm=transform_to_log_scale(data.hASP$TP.ASP-mean(na.omit(data.hASP$TP.ASP))))
data.hASP =dplyr::rename(data.hASP, TP = TP.ASP)

data.hPRO = data[,c("years","Location.2","TP.PRO")]
data.hPRO =cbind(data.hPRO, AA=rep('Pro',length(data.hPRO$TP.PRO)), TP.norm=transform_to_log_scale(data.hPRO$TP.PRO-mean(na.omit(data.hPRO$TP.PRO))))
data.hPRO =dplyr::rename(data.hPRO, TP = TP.PRO)

data.hVAL = data[,c("years","Location.2","TP.VAL")]
data.hVAL =cbind(data.hVAL, AA=rep('Val',length(data.hVAL$TP.VAL)), TP.norm=transform_to_log_scale(data.hVAL$TP.VAL-mean(na.omit(data.hVAL$TP.VAL))))
data.hVAL =dplyr::rename(data.hVAL, TP = TP.VAL)

data.hALA = data[,c("years","Location.2","TP.ALA")]
data.hALA =cbind(data.hALA, AA=rep('Ala',length(data.hALA$TP.ALA)), TP.norm=transform_to_log_scale(data.hALA$TP.ALA-mean(na.omit(data.hALA$TP.ALA))))
data.hALA =dplyr::rename(data.hALA, TP = TP.ALA)

########################     Hierarchical Environmental Models       ############################

data.hier <- rbind(data.hGLU, data.hASP)
data.hier <- rbind(data.hier, data.hPRO)
data.hier <- rbind(data.hier, data.hVAL)
data.hier <- rbind(data.hier, data.hALA)
data.hier<- subset(data.hier, Location.2=="Coastal"|Location.2=="Inland")
data.hier <- cbind(data.hier, Year=data.hier$years+1)

dfa <- read.csv("Data/Compiled/DFA.csv")
dfa = dfa[,c("Year","Up.DFA1","Climate.DFA1.x", "SST.DFA1", "Dis.DFA1")]
data2<-merge(dfa,data.hier, by='Year', all.x=TRUE)


Mod.Upwelling <- lmer(TP~Up.DFA1+Location.2+(1|AA), data=data2)
summary(Mod.Upwelling)
AIC(Mod.Upwelling)

Mod <- lmer(TP~Location.2+(1|AA), data=data2)
summary(Mod)
AIC(Mod)

Mod.Clim <- lmer(TP~Climate.DFA1.x+Location.2+(1|AA), data=data2)
summary(Mod.Clim)
AIC(Mod.Clim)


ModelSelection.WA <- function(dataframe,n,y) {
  
  aic.output <- rbind(AIC(lmer(y~Location.2+(1|AA), data=dataframe)), 
                      AIC(lmer(y~Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~SST.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Up.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Dis.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(y~Climate.DFA1.x+SST.DFA1+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(y~Climate.DFA1.x+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(y~Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(y~Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(y~Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(y~Up.DFA1+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#7
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#8
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#9
                      AIC(lmer(y~Climate.DFA1.x+Up.DFA1+Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#10
                      AIC(lmer(y~Climate.DFA1.x+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe))#11
  )
  
  model.summary <- rbind(summary(lmer(y~Location.2+(1|AA), data=dataframe)), 
                         summary(lmer(y~Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(y~SST.DFA1+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(y~Up.DFA1+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(y~Dis.DFA1+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(y~Climate.DFA1.x+SST.DFA1+Location.2+(1|AA), data=dataframe)),#1
                         summary(lmer(y~Climate.DFA1.x+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#2
                         summary(lmer(y~Climate.DFA1.x+Up.DFA1+Location.2+(1|AA), data=dataframe)),#3
                         summary(lmer(y~Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#4
                         summary(lmer(y~Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#5
                         summary(lmer(y~Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#6
                         summary(lmer(y~Up.DFA1+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#7
                         summary(lmer(y~Climate.DFA1.x+Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#8
                         summary(lmer(y~Climate.DFA1.x+Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#9
                         summary(lmer(y~Climate.DFA1.x+Up.DFA1+Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#10
                         summary(lmer(y~Climate.DFA1.x+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe))#11
                         
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Climate", "SST", "Up", "Dis","1. Clim, SST", "2. Clim, Dis", "3. Clim, Up", "4. SST, Dis",
                   "5. SST, Up", "6. Up, Dis", "7. SST, Dis, Up", "8. Clim, Up, Dis", "9. Clim, Up, SST", "10. Clim, Up, SST, Dis", "11. Clim Dis SST")
  
  row.names(aic.output) <- model.names
  delaic <- aic.output-min(aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  aic.output <- cbind(aic.output, delaic, aic.weight)
  aic.output<- cbind(aic.output, model.names)
  colnames(aic.output)<- c("AICc", "delAICc", "AICc Weight", "model.names")
  return(aic.output)
}


n<- 17
model.selectionENV <- ModelSelection.WA(data2, n, data2$TP.norm)
model.selectionENV <-transform(model.selectionENV, AICc = as.numeric(AICc))

########################     Hierarchical Prey Models       ############################

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
data3<-merge(prey,data.hier, by='Year', all.x=TRUE)
data3 <- subset(data3, Year>=1973&Year<=2008)


Mod.loc <- lmer(TP~Location.2+(1|AA), data=data3)
summary(Mod.loc)
AIC(Mod.loc)

Mod <- lmer(TP~Location.2+Herring.Biomass+(1|AA), data=data3)
summary(Mod)
AIC(Mod)

Mod <- lmer(TP~Location.2+Chinook+(1|AA), data=data3)
summary(Mod)
AIC(Mod)

Mod <- lmer(TP~Location.2+HatcherySmolts+(1|AA), data=data3)
summary(Mod)
AIC(Mod)

Mod <- lmer(TP~Location.2+HakeBiomass+(1|AA), data=data3)
summary(Mod)
AIC(Mod)
