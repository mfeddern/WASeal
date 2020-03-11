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


data.hGLU = data[,c("years","Location.2","TP.GLU", "Sample.ID")]
data.hGLU =cbind(data.hGLU, AA=rep('Glu',length(data.hGLU$TP.GLU)), TP.norm=data.hGLU$TP.GLU-mean(na.omit(data.hGLU$TP.GLU)))
data.hGLU =dplyr::rename(data.hGLU, TP = TP.GLU)

data.hASP = data[,c("years","Location.2","TP.ASP", "Sample.ID")]
data.hASP =cbind(data.hASP, AA=rep('Asp',length(data.hASP$TP.ASP)), TP.norm=data.hASP$TP.ASP-mean(na.omit(data.hASP$TP.ASP)))
data.hASP =dplyr::rename(data.hASP, TP = TP.ASP)

data.hPRO = data[,c("years","Location.2","TP.PRO", "Sample.ID")]
data.hPRO =cbind(data.hPRO, AA=rep('Pro',length(data.hPRO$TP.PRO)), TP.norm=data.hPRO$TP.PRO-mean(na.omit(data.hPRO$TP.PRO)))
data.hPRO =dplyr::rename(data.hPRO, TP = TP.PRO)

data.hVAL = data[,c("years","Location.2","TP.VAL", "Sample.ID")]
data.hVAL =cbind(data.hVAL, AA=rep('Val',length(data.hVAL$TP.VAL)), TP.norm=data.hVAL$TP.VAL-mean(na.omit(data.hVAL$TP.VAL)))
data.hVAL =dplyr::rename(data.hVAL, TP = TP.VAL)

data.hALA = data[,c("years","Location.2","TP.ALA", "Sample.ID")]
data.hALA =cbind(data.hALA, AA=rep('Ala',length(data.hALA$TP.ALA)), TP.norm=data.hALA$TP.ALA-mean(na.omit(data.hALA$TP.ALA)))
data.hALA =dplyr::rename(data.hALA, TP = TP.ALA)

data.hPHE = data[,c("years","PHE.mean", "Sample.ID")]
data.hPHE =cbind(data.hPHE,PHE.norm=data.hPHE$PHE.mean-mean(na.omit(data.hPHE$PHE.mean)))

data.d13C = data[,c("years","d13C.s", "Sample.ID")]
data.d13C =cbind(data.d13C,d13C.norm=data.d13C$d13C.s-mean(na.omit(data.d13C$d13C.s)))

########################     Hierarchical Climate Models       ############################

data.hier <- rbind(data.hGLU, data.hASP)
data.hier <- rbind(data.hier, data.hPRO)
data.hier <- rbind(data.hier, data.hVAL)
data.hier <- rbind(data.hier, data.hALA)
data.hier<- subset(data.hier, Location.2=="Coastal"|Location.2=="Inland")
data.hier <- cbind(data.hier, Year=data.hier$years+1) # adding 1 year lag

dfa <- read.csv("Data/Compiled/DFA.csv")
dfa = dfa[,c("Year","Up.DFA1","Climate.DFA1.x", "SST.DFA1", "Dis.DFA1")]
data2<-merge(dfa,data.hier, by='Year', all.x=TRUE)
data2<-merge(data2,data.hPHE, by='Sample.ID', all.x=TRUE)
data2<-merge(data2,data.d13C, by='Sample.ID', all.x=TRUE)


ModelSelection.WA <- function(dataframe,n) {
  
  aic.output <- rbind(AIC(lmer(TP.norm~Location.2+(1|AA), data=dataframe)), 
                      AIC(lmer(TP.norm~Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~SST.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~Up.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~Dis.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~Climate.DFA1.x+SST.DFA1+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(TP.norm~Climate.DFA1.x+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(TP.norm~Climate.DFA1.x+Up.DFA1+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(TP.norm~Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(TP.norm~Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(TP.norm~Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(TP.norm~Up.DFA1+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#7
                      AIC(lmer(TP.norm~Climate.DFA1.x+Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#8
                      AIC(lmer(TP.norm~Climate.DFA1.x+Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#9
                      AIC(lmer(TP.norm~Climate.DFA1.x+Up.DFA1+Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#10
                      AIC(lmer(TP.norm~Climate.DFA1.x+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe))
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
model.selectionENV <- ModelSelection.WA(data2, n)
model.selectionENV <-transform(model.selectionENV, AICc = as.numeric(AICc))

########################     Hierarchical Prey Models       ############################

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
data.merge <- list(prey,data.hier, dfa)
data3 <- Reduce(function(x,y) merge(x,y, by='Year', all.x=TRUE), data.merge)
data3<-merge(data3,data.hPHE, by='Sample.ID', all.x=TRUE)
data3<-merge(data3,data.d13C, by='Sample.ID', all.x=TRUE)
data3 <- subset(data3, Year>=1973&Year<=2008)


ModelSelection.WA2 <- function(dataframe,n) {
  
  aic.output <- rbind(AIC(lmer(TP.norm~Location.2+(1|AA), data=dataframe)), 
                      AIC(lmer(TP.norm~Herring.Biomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~Chinook+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~HatcherySmolts+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~HakeBiomass+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~Herring.Biomass+Chinook+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(TP.norm~Herring.Biomass+HakeBiomass+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(TP.norm~Herring.Biomass+HatcherySmolts+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(TP.norm~HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(TP.norm~HatcherySmolts+Chinook+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(TP.norm~HatcherySmolts+HakeBiomass+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(TP.norm~HatcherySmolts+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)),#7
                      AIC(lmer(TP.norm~Herring.Biomass+HatcherySmolts+HakeBiomass+Location.2+(1|AA), data=dataframe)),#8
                      AIC(lmer(TP.norm~Herring.Biomass+HatcherySmolts+Chinook+Location.2+(1|AA), data=dataframe)),#9
                      AIC(lmer(TP.norm~Herring.Biomass+HatcherySmolts+HakeBiomass+Chinook+Location.2+(1|AA), data=dataframe)),#10
                      AIC(lmer(TP.norm~Herring.Biomass+Chinook+HakeBiomass+Location.2+(1|AA), data=dataframe)) #11
                      
            
                      
  )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
                   "5. Chinook, Hatch", "6. Hatch, Hake", "7. Chinook, Hake, Hatch", "8. Herring, Hatch, Hake", "9. Herring, Hatch, Chinook", "10. Herring, 
                   Hatch, Chinook, Hake", "11. Herring Hake Chinook")
  
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

model.selectionPREY <- ModelSelection.WA2(data3, n)
model.selectionPREY <-transform(model.selectionPREY, AICc = as.numeric(AICc))




########################     Hierarchical Nutrient Models       ############################


ModelSelection.WA3 <- function(dataframe,n) {
  
  aic.output <- rbind(AIC(lmer(TP.norm~Location.2+(1|AA), data=dataframe)), 
                      AIC(lmer(TP.norm~PHE.norm+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~d13C.norm+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP.norm~PHE.norm+d13C.norm+Location.2+(1|AA), data=dataframe))
                      
              
                       )
  
  names <- seq(1,n,1)
  model.names <- c("Location", "Phe", "13C", "PHE, 13C")
  
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
model.selectionNUTRIENT <- ModelSelection.WA3(data2, n)
model.selectionPREY <-transform(model.selectionPREY, AICc = as.numeric(AICc))

#AIC(lmer(TP.norm~PHE.mean+(1|AA), data=dataframe)),#10




########################     Lasso Regression       ############################


WA.lasso <- data3
cor(WA.lasso)
WA.lasso <- subset(WA.lasso, AA="Glu")
Wa.lasso <- WA.lasso[complete.cases(WA.lasso), ]
WA.lasso <- na.omit(WA.lasso)

y<- Wa.lasso$TP.norm
Wa.lasso <- Wa.lasso %>%
  select(-TP) %>%
  select(-Sample.ID) %>%
  select(-Year) %>%
  select(-X) %>%
  select(-years.x) %>%
  select(-years.y) %>%
  select(-TP.norm) %>%
  select(-AA)%>%
  select(-PHE.mean)%>%
  select(-d13C.s)%>%
  select(-years)%>%
  #select(-PHE.norm)%>%
  select(-Climate.DFA1.x)%>%
  select(-HarborSeal)
  #select(-Climate.DFA1.x)%>%
  #select(-Up.DFA1)%>%
  #select(-SST.DFA1)%>%
  #select(-Dis.DFA1)


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
numLevels[numLevels==0] <- 1
#apply(X[, 'Location'], factor)

#X[8][X[8] == 2] <- 0
#X[8][X[8] == 3] <- 1

# make the categorical variables take integer values starting from 0
X[, !i_num] <- apply(X[, !i_num], 2, function(col) as.integer(as.factor(col)) - 1)
sapply(X, class)






library(glinternet)
set.seed(2001)

head(X)

pairs <- matrix(c(8,1, 8,2, 8,3, 8,4,
                  8,5, 8,6, 8,7,  8,9, 
                  8,10,  8,11,  8,12, 8,13, 8,14), ncol=13, byrow=TRUE)

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



