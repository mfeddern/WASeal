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

data.hGLU = data[,c("years","Location.2","TP.GLU")]
data.hGLU =cbind(data.hGLU, AA=rep('Glu',length(data.hGLU$TP.GLU)))
data.hGLU =dplyr::rename(data.hGLU, TP = TP.GLU)

data.hASP = data[,c("years","Location.2","TP.ASP")]
data.hASP =cbind(data.hASP, AA=rep('Asp',length(data.hASP$TP.ASP)))
data.hASP =dplyr::rename(data.hASP, TP = TP.ASP)

data.hPRO = data[,c("years","Location.2","TP.PRO")]
data.hPRO =cbind(data.hPRO, AA=rep('Pro',length(data.hPRO$TP.PRO)))
data.hPRO =dplyr::rename(data.hPRO, TP = TP.PRO)

data.hVAL = data[,c("years","Location.2","TP.VAL")]
data.hVAL =cbind(data.hVAL, AA=rep('Val',length(data.hVAL$TP.VAL)))
data.hVAL =dplyr::rename(data.hVAL, TP = TP.VAL)

data.hALA = data[,c("years","Location.2","TP.ALA")]
data.hALA =cbind(data.hALA, AA=rep('Ala',length(data.hALA$TP.ALA)))
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


ModelSelection.WA <- function(dataframe,n) {
  
  aic.output <- rbind(AIC(lmer(TP~Location.2+(1|AA), data=dataframe)), 
                      AIC(lmer(TP~Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP~SST.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP~Up.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP~Dis.DFA1+Location.2+(1|AA), data=dataframe)),
                      AIC(lmer(TP~Climate.DFA1.x+SST.DFA1+Location.2+(1|AA), data=dataframe)),#1
                      AIC(lmer(TP~Climate.DFA1.x+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#2
                      AIC(lmer(TP~Climate.DFA1.x+Up.DFA1+Location.2+(1|AA), data=dataframe)),#3
                      AIC(lmer(TP~Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#4
                      AIC(lmer(TP~Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#5
                      AIC(lmer(TP~Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#6
                      AIC(lmer(TP~Up.DFA1+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#7
                      AIC(lmer(TP~Climate.DFA1.x+Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#8
                      AIC(lmer(TP~Climate.DFA1.x+Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#9
                      AIC(lmer(TP~Climate.DFA1.x+Up.DFA1+Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#10
                      AIC(lmer(TP~Climate.DFA1.x+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe))#11
  )
  
  model.summary <- rbind(summary(lmer(TP~Location.2+(1|AA), data=dataframe)), 
                         summary(lmer(TP~Climate.DFA1.x+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(TP~SST.DFA1+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(TP~Up.DFA1+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(TP~Dis.DFA1+Location.2+(1|AA), data=dataframe)),
                         summary(lmer(TP~Climate.DFA1.x+SST.DFA1+Location.2+(1|AA), data=dataframe)),#1
                         summary(lmer(TP~Climate.DFA1.x+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#2
                         summary(lmer(TP~Climate.DFA1.x+Up.DFA1+Location.2+(1|AA), data=dataframe)),#3
                         summary(lmer(TP~Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#4
                         summary(lmer(TP~Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#5
                         summary(lmer(TP~Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#6
                         summary(lmer(TP~Up.DFA1+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#7
                         summary(lmer(TP~Climate.DFA1.x+Up.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe)),#8
                         summary(lmer(TP~Climate.DFA1.x+Up.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#9
                         summary(lmer(TP~Climate.DFA1.x+Up.DFA1+Dis.DFA1+SST.DFA1+Location.2+(1|AA), data=dataframe)),#10
                         summary(lmer(TP~Climate.DFA1.x+SST.DFA1+Dis.DFA1+Location.2+(1|AA), data=dataframe))#11
                         
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


n<- 16
model.selectionENV <- ModelSelection.WA(data2, n)
model.selectionENV <-transform(model.selectionENV, AICc = as.numeric(AICc))

#################################################################################################
########################     Normalizing and transforming Seal data ############################
#################################################################################################


Inland <- subset(data, Location.2=="Inland")
Coastal <- subset(data, Location.2=="Coastal")


#data <- cbind(data[,'Year'],data[,'TP.GLU'], data[,'PHE.mean'])

Coastal <- cbind(Coastal[,'Year'],Coastal[,'TP.GLU'], Coastal[,'PHE.mean'])
Inland <- cbind(Inland[,'Year'],Inland[,'TP.GLU'], Inland[,'PHE.mean'])


aggdataC <-aggregate(Coastal, by=list(Coastal[,1]), 
                     FUN=mean, na.rm=TRUE)
aggdataI <-aggregate(Inland, by=list(Inland[,1]), 
                     FUN=mean, na.rm=TRUE)


TP.C <- aggdataC[,2:4]
TP.I <- aggdataI[,2:4]

colnames(TP)<- c("Year", "TPC", "PHEC", "TPI", "PHEI")

logdatI <- log(TP.I[,2:3])
logdatC <- log(TP.C[,2:3])
logdatC <- logdatC[-7,]

y_barTPC <- mean(logdatC[,1])
y_barPHEC <- mean(logdatC[,2])
y_barTPI <- mean(logdatI[,1])
y_barPHEI <- mean(logdatI[,2])
c.year<-TP.C[,1]
c.year<- c.year[-7]

log.normC <- cbind(c.year,c(logdatC[,1] - y_barTPC),c(logdatC[,2] - y_barPHEC))
log.normI <- cbind(TP.I[,1],c(logdatI[,1] - y_barTPI),c(logdatI[,2] - y_barPHEI))

#logdatC <- cbind(c.year,logdatC)
#logdatI <- cbind(TP.I[,1],logdatI)

#log.normC <- cbind(c.year,c(logdatC[,1] - y_barTPC),c(logdatC[,2] - y_barPHEC))
#log.normI <- cbind(TP.I[,1],c(logdatI[,1] - y_barTPI),c(logdatI[,2] - y_barPHEI))

logdatC<- log.normC
logdatI<- log.normI


colnames(logdatC) <- c("Year", "TPC", "PHEC")
colnames(logdatI) <- c("Year", "TPI", "PHEI")


#################################################################################################
########################                All WA Salmon                 ############################
#################################################################################################
WAsalm <- read.csv("Salmon/Fraser.csv")
head(WAsalm)
WAsalm.tot <- WAsalm$Total
log.was <- log(WAsalm.tot)  
y_bar <- mean(log.was) 
was.norm <- log.was-y_bar
WAsalm <- cbind(WAsalm$Year,was.norm)
colnames(WAsalm)<- c("Year", "SalmWA")
plot(WAsalm[,1], WAsalm[,2], type='l')


dat <-merge(WAsalm, logdatC, by='Year', all=TRUE)
dat <-merge(dat, logdatI, by='Year', all=TRUE)

#################################################################################################
########################                Fitting the MARSS Model TP                 ############################
#################################################################################################



plot(dat[,1], dat[,3], ylim = c(-1,1), col='green3', pch=16, ylab='log(Trophic Position)', xlab='Year')
points(dat[,1], dat[,5], col='blue', pch=16)


plot(dat[,'Year'], dat[,'PHEC'], col='green3', pch=16, ylab='log(PHE d15N)', xlab='Year', ylim=c(-1,1))
points(dat[,'Year'], dat[,'PHEI'], col='blue', pch=16)

plot(dat[,1], dat[,3], ylim = c(-1,1), col='green3', pch=16, ylab='Trophic Position', xlab='Year')
points(dat[,1], dat[,5], col='blue', pch=16)


seal.dat <- cbind(dat[,'TPC'], dat[,'TPI'])
years = dat[,1]
row.names(seal.dat) <- years
dat.1 = t(seal.dat)  #transpose to have years across columns
n = nrow(dat.1) - 1

Z=matrix(c(1,1)) ###writing z matrix for single population
mod.list.0 = list(B = matrix(1), U = matrix("u"), Q ="diagonal and unequal", 
                  Z = Z, R = "diagonal and equal", A="scaling",
                  x0 = matrix("mu"), tinitx = 0)
fit.1 = MARSS(dat.1, model = mod.list.0) #-79.15746 




Z=matrix(c(1,0,
           0,1),nrow=2, byrow = TRUE) ###writing z matrix for two populations
B <- matrix(c(1,0,
              0,1), nrow=2)
U <- matrix(c("uC", "uI"))
x0 <- matrix(c("C", "I"))

mod.list.0 = list(B = B, U = U, Q = "diagonal and equal", 
                  Z = Z, R = "diagonal and equal", A="scaling",
                  x0 = x0, tinitx = 0)

fit.2 = MARSS(dat.1, model = mod.list.0) #AICc -73.18695


mod.list.0 = list(B = B, U = "equal", Q = "diagonal and equal", 
                  Z = Z, R = "diagonal and equal", A="scaling",
                  x0 = x0, tinitx = 0)

fit.8 = MARSS(dat.1, model = mod.list.0) #AICc -73.18695

mod.list.0 = list(B = B, U = U, Q = "diagonal and unequal", 
                  Z = Z, R = "diagonal and equal", A="scaling",
                  x0 = x0, tinitx = 0)

fit.3 = MARSS(dat.1, model = mod.list.0) #AICc -73.3307****convergence issues

zsalm <- zscore(dat[,'SalmWA'])
salm.dat <- t(dat[,'SalmWA'])
fit.7 = MARSS(dat.2, model = mod.list.0, covariates = salm.dat) #AICc -34.51 ***convergence issues

#################################################################################################
########################                Plotting Residuals                ############################
#################################################################################################



resids <- residuals(fit.2)$model.residuals
#resids[is.na(dat)] <- NA



plot(resids[1,], pch=16, ylab="Model Residuals")
abline(h=0)
plot(resids[2,], pch=16, ylab="Model Residuals")
abline(h=0)
hist(resids[1,])
hist(resids[2,])

acf(resids[1,], na.action=na.pass)
acf(resids[2,], na.action=na.pass)

auto.arima(resids[1,])
auto.arima(resids[2,])


#################################################################################################
########################                Fitting the MARSS Model PHE              ############################
#################################################################################################




seal.dat <- cbind(dat[,'PHEC'], dat[,'PHEI'])
years = seal.dat[,1]
dat.2 = seal.dat[, !(colnames(seal.dat) %in% c("Year"))]
dat.2 = t(seal.dat)  #transpose to have years across columns
colnames(dat.2) = years
n = nrow(dat.2) - 1

Z=matrix(c(1,1)) ###writing z matrix for single population
mod.list.0 = list(B = matrix(1), U = matrix("u"), Q = matrix("q"), 
                  Z = Z, R = "diagonal and equal", A="scaling",
                  x0 = matrix("mu"), tinitx = 0)
fit.4 = MARSS(dat.2, model = mod.list.0) #-39.34 AICc




Z=matrix(c(1,0,
           0,1),nrow=2, byrow = TRUE) ###writing z matrix for two populations
B <- matrix(c(1,0,
              0,1), nrow=2)
U <- matrix(c('uC',"uI"))
x0 <- matrix(c("C", "I"))

mod.list.0 = list(B = B, U = U, Q = "diagonal and equal", 
                  Z = Z, R = "diagonal and equal", A="scaling",
                  x0 = x0, tinitx = 0)

fit.5 = MARSS(dat.2, model = mod.list.0) #AICc -37.58333

mod.list.0 = list(B = B, U = U, Q = "diagonal and unequal", 
                  Z = Z, R = "diagonal and equal", A="scaling",
                  x0 = x0, tinitx = 0)

fit.6 = MARSS(dat.2, model = mod.list.0) #AICc -34.51 ***convergence issues


#################################################################################################
########################                Plotting Residuals                ############################
#################################################################################################



resids <- residuals(fit.5)$model.residuals
#resids[is.na(dat)] <- NA



plot(resids[1,], pch=16, ylab="Model Residuals")
abline(h=0)
plot(resids[2,], pch=16, ylab="Model Residuals")
abline(h=0)
hist(resids[1,])
hist(resids[2,])

acf(resids[1,], na.action=na.pass)
acf(resids[2,], na.action=na.pass)

auto.arima(resids[1,])
auto.arima(resids[2,])
#################################################################################################
########################              DFA                ############################
#################################################################################################

yr_frst <- 1975
yr_last <- 2017
dat <- dat[dat[,1]>=yr_frst & dat[,1]<=yr_last,]

plot(dat[,1], dat[,2],type = 'l')
lines(dat[,1], dat[,3], col='red')
lines(dat[,1], dat[,4], col='blue')
lines(dat[,1], dat[,5], col='green')
lines(dat[,1], dat[,6], col='purple')

dat2<-dat[,3:6]
row.names(dat2) <- dat[,1]

dat.t <- t(dat2)
dat.t <- rbind(dat.t[2,], dat.t[4,])
row.names(dat.t) <- c("TPC", "TPI")
N_ts <- dim(dat.t )[1]
## get length of time series
TT <- dim(TP.t)[2] 

aa <- "zero"
DD <- "zero"  # matrix(0,mm,1)
dd <- "zero"  # matrix(0,1,wk_last)
RR <-"diagonal and equal"
#RR<- matrix(0, ncol = n, nrow = n)


type <- rownames(dat.t)
clr <- c("brown", "blue", "darkgreen", "darkred", "purple")
cnt <- 1
ncomp <- 1
##DFA with 1 trend
n<-2
mm<- 1
Z_vals <- list("z11",
               "z21")

ZZ <- matrix(Z_vals, nrow=n, ncol=mm, byrow=TRUE)
ZZ

BB <- "identity"  # diag(mm)
uu <- "zero"  # matrix(0,mm,1)
CC <- "zero"  # matrix(0,mm,1)
cc <- "zero"  # matrix(0,1,wk_last)
QQ <- "identity"  # diag(mm)


mod_list <- list(Z=ZZ, A=aa, D=DD, d=dd, R=RR,
                 B=BB, U=uu, C=CC, c=cc, Q=QQ)
init_list <- list(x0 = matrix(rep(0, mm), mm, 1))
con_list <- list(maxit = 3000, allow.degen = TRUE)

#dfa_1 <- MARSS(y = dat.t, model = mod_list, inits = init_list, control = con_list)
#dfa_2 <- MARSS(y = dat.t, model = mod_list, inits = init_list, control = con_list)
dfa_3 <- MARSS(y = dat.t, model = mod_list, inits = init_list, control = con_list)


ylbl <- c("TPC", "TPI")
## get the estimated ZZ
Z_est <- coef(dfa_3, type="matrix")$Z
## Code for plotting one trend 

Z_est <- coef(dfa_3, type="matrix")$Z
H_inv <- 1
Z_rot = Z_est %*% H_inv   
proc_rot = solve(H_inv) %*% dfa_3$states
ylbl <- phytoplankton
w_ts <- seq(dim(dat.t)[2])
layout(matrix(c(1,2),mm,2),widths=c(2,1))
par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
## plot the processes
for(i in 1:mm) {
  ylm <- c(-1,1)*max(abs(proc_rot[i,]))
  ## set up plot area
  plot(w_ts,proc_rot[i,], type="n", bty="L",
       ylim=ylm, xlab="", ylab="", xaxt="n")
  ## draw zero-line
  abline(h=0, col="gray")
  ## plot trend line
  lines(w_ts,proc_rot[i,], lwd=2)
  lines(w_ts,proc_rot[i,], lwd=2)
  ## add panel labels
  mtext(paste("State",i), side=3, line=0.5)
  axis(1,12*(0:dim(dat.t)[2])+1,yr_frst+0:dim(dat.t)[2])
}
## plot the loadings
minZ <- 0
ylm <- c(-1,1)*max(abs(Z_rot))
for(i in 1:mm) {
  plot(c(1:N_ts)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
       lwd=2, xlab="", ylab="", xaxt="n", ylim=ylm, xlim=c(0.5,N_ts+0.5), col=clr)
  for(j in 1:N_ts) {
    if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
    if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
    abline(h=0, lwd=1.5, col="gray")
  } 
  mtext(paste("Factor loadings on state",i),side=3,line=0.5)
  
}





Z_est <- coef(dfa_3, type="matrix")$Z
H_inv <- varimax(Z_est)$rotmat
Z_rot = Z_est %*% H_inv   
proc_rot = solve(H_inv) %*% dfa_3$states


ylbl <- c("PDO","TP", "PHE")

w_ts <- seq(dim(dat.t)[2])
layout(matrix(c(1, 2, 3, 4, 5, 6), mm, 2), widths = c(2, 1))
## par(mfcol=c(mm,2), mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0))
par(mai = c(0.5, 0.5, 0.5, 0.1), omi = c(0, 0, 0, 0))
## plot the processes
for (i in 1:mm) {
  ylm <- c(-1, 1) * max(abs(proc_rot[i, ]))
  ## set up plot area
  plot(w_ts, proc_rot[i, ], type = "n", bty = "L", ylim = ylm, 
       xlab = "", ylab = "", xaxt = "n")
  ## draw zero-line
  abline(h = 0, col = "gray")
  ## plot trend line
  lines(w_ts, proc_rot[i, ], lwd = 2)
  lines(w_ts, proc_rot[i, ], lwd = 2)
  ## add panel labels
  mtext(paste("State", i), side = 3, line = 0.5)
  axis(1, 0:dim(dat.t)[2] + 1, yr_frst + 0:dim(dat.t)[2])
}







## plot the loadings
minZ <- 0
ylm <- c(-1, 1) * max(abs(Z_rot))
for (i in 1:mm) {
  plot(c(1:n)[abs(Z_rot[, i]) > minZ], as.vector(Z_rot[abs(Z_rot[, 
                                                                 i]) > minZ, i]), type = "h", lwd = 2, xlab = "", ylab = "", 
       xaxt = "n", ylim = ylm, xlim = c(0.5, N_ts + 0.5), col = clr)
  for (j in 1:n) {
    if (Z_rot[j, i] > minZ) {
      text(j, -0.03, ylbl[j], srt = 90, adj = 1, cex = 1.2, 
           col = clr[j])
    }
    if (Z_rot[j, i] < -minZ) {
      text(j, 0.03, ylbl[j], srt = 90, adj = 0, cex = 1.2, 
           col = clr[j])
    }
    abline(h = 0, lwd = 1.5, col = "gray")
  }
  mtext(paste("Factor loadings on state", i), side = 3, line = 0.5)
}




################################################################################################################################
##################################################           Running Hier Models  ##################################################
################################################################################################################################
hier<-read.csv("Hierarchical.csv")

h.year <- lmer(d15N~years+(1|AA), data=hier)
summary(h.year)
AIC(h.year) #3483.93

h.year2 <- lmer(d15N~years+(1+years|AA), data=hier)
summary(h.year2)
AIC(h.year2) #3483.93

h.year3 <- lm(PHE.mean~years, data=dat1)
summary(h.year3)
AIC(h.year3) #979.1555


h.location <- lmer(d15N~Location.2+(1|AA), data=hier)
summary(h.location)
AIC(h.location) #3435.396

h.location2 <- lmer(d15N~Location.2+(1+Location.2|AA), data=hier)
summary(h.location2)
AIC(h.location2) #3435.39

h.location3 <- lm(PHE.mean~Location.2, data=dat1)
summary(h.location3)
AIC(h.location3) #924.1471


h.ly <- lmer(d15N~Location.2:years+(1|AA), data=hier)
summary(h.ly)
AIC(h.ly) #3503.65

h.ly <- lmer(d15N~years+Location.2+(1|AA), data=hier)
summary(h.ly)
AIC(h.ly) #3503.65

h.ly2 <- lmer(d15N~years+Location.2+(1+years|AA), data=hier)
summary(h.ly2)
AIC(h.ly) #3503.65


h.ly3 <- lm(PHE.mean~years, data=dat1)
summary(h.ly3)
AIC(h.ly3) #926.143

h.ly3 <- lm(PHE.mean~Location.2, data=dat1)
summary(h.ly3)
boxplot(h.ly3)
AIC(h.ly3) #926.143
summary(aov(h.ly3))


subset(dat1,PHE.mean >= 20)
boxplot(PHE.mean~Location.2, frame.plot=FALSE,data=dat1,
        xlab="Location", pch=16,  col=c("#333399","#999966", "goldenrod", "#006699", "#993399"), ylim=c(0,25), xaxt="n")
mtext(text=expression(paste(delta^{15}, "N (‰) Phenylalanine")), 2, 2.5, bty='n')
text(1, -1, labels="Eastern", pos=1, xpd=TRUE, cex=0.87)
text(1, -2, labels="Bering Sea", pos=1, xpd=TRUE, cex=0.87)
text(2, -1, labels="Coastal", pos=1, xpd=TRUE, cex=0.87)
text(2, -2, labels="WA", pos=1, xpd=TRUE, cex=0.87)
text(3, -1, labels="Salish Sea", pos=1, xpd=TRUE, cex=0.87)
text(4, -1, labels="Southcentral", pos=1, xpd=TRUE, cex=0.87)
text(4, -2, labels="AK", pos=1, xpd=TRUE, cex=0.87)
text(5, -1, labels="Southeast", pos=1, xpd=TRUE, cex=0.87)
text(5, -2, labels="AK", pos=1, xpd=TRUE, cex=0.87)

################################################################################################################################
##################################################           Running GAMS   AK   ##################################################
################################################################################################################################

datGOA <- subset(dat1, Location.2=="SE"|Location.2=="SC")
fit.1 <- lm(PHE.mean~TUMI.45.lag+Location.2, data=datGOA)
summary(fit.1)
AIC(fit.1)

################################################################################################################################
##################################################          Washington      #########################################
################################################################################################################################

datWA <- subset(dat1, Location.2=="Coastal"|Location.2=="Inland")
fit.1 <- lm(PHE.mean~MEI.YEAR.lag+Location.2, data=datWA)
summary(fit.1)
AIC(fit.1)


fit.1 <- lm(PHE.mean~TUMI.45.lag+Location.2, data=datWA)
summary(fit.1)
AIC(fit.1)

fit.1 <- lm(PHE.mean~TUMI.45.lag+MEI.YEAR.lag+Location.2, data=datWA)
summary(fit.1)
AIC(fit.1)


fit.1 <- lm(PHE.mean~MEI.YEAR.lag+Location.2, data=datWA)
summary(fit.1)
AIC(fit.1)#434.1876

fit.1 <- lm(PHE.mean~MEI.YEAR.lag*Location.2, data=datWA)
summary(fit.1)
AIC(fit.1)#435.3356

fit.1 <- lm(PHE.mean~MEI.YEAR.lag+TUMI.45.lag+Location.2, data=datWA)
summary(fit.1)
AIC(fit.1)#431.5327


plot(datWA$PDO.lag, datWA$PHE.mean, ylim=c(4, 22))
abline(11.7888,0.4065)

fit.1 <- gam(PHE.mean~s(MEI.YEAR.lag), data=dat2)
summary(fit.1)


datSE <- subset(dat1, Location.2=="SE")
fit.1 <- lm(PHE.mean~MEI.YEAR.lag, data=datSC)
summary(fit.1)
AIC(fit.1)
plot(datSC$MEI.YEAR.lag, datSC$PHE.mean, ylim=c(4, 22))


datSC <- subset(dat1, Location.2=="SC")
fit.1 <- lm(PHE.mean~SSH.YEAR.lag, data=datSC)
summary(fit.1)
AIC(fit.1)
plot(datSC$MEI.YEAR.lag, datSC$PHE.mean, ylim=c(4, 22))

datBB <- subset(dat1, Location.2=="BB")
fit.1 <- lm(PHE.mean~SSH.YEAR.lag, data=datBB)
summary(fit.1)
AIC(fit.1)
plot(datSC$MEI.YEAR.lag, datSC$PHE.mean, ylim=c(4, 22))

















fit.1 <- gam(PHE.mean~Location.2, data=dat1)
summary(fit.1)
AIC(fit.1)

fit.1 <- gam(PHE.mean~years, data=dat1)
summary(fit.1)
AIC(fit.1)

fit.1 <- gam(PHE.mean~years:Location.2, data=dat1)
summary(fit.1)
AIC(fit.1)

fit.1 <- gam(PHE.mean~Location.2, data=dat1)
summary(fit.1)
AIC(fit.1) #924.1471

fit.1 <- gam(PHE.mean~MEI.YEAR.lag*Location.2, data=dat1)
summary(fit.1)
AIC(fit.1) #434.1876; 869.3463; 872.4981

fit.1 <- gam(PHE.mean~MEI.YEAR.lag*Location.2+PDO.lag*Location.2, data=dat1)
summary(fit.1)
AIC(fit.1) #435.2366; 871.2905; 

fit.1 <- gam(PHE.mean~MEI.YEAR.lag*Location.2+TUMI.45.lag*Location.2, data=dat1)
summary(fit.1)
AIC(fit.1) #431.5327; 811.0825*****BEST 819.1489

fit.1 <- gam(PHE.mean~MEI.YEAR.lag+PDO.lag+TUMI.45.lag+Location.2, data=dat1)
summary(fit.1)
AIC(fit.1) #433.3691; 813.0313





fit.1 <- gam(PHE.mean~PDO.lag+Location.2, data=dat1)
summary(fit.1)
AIC(fit.1) #480.4504; 923.9078

fit.1 <- gam(PHE.mean~PDO.lag*Location.2+TUMI.45.lag*Location.2, data=dat1)
summary(fit.1)
AIC(fit.1) #432.9861; 811.8282*****STILL GOOD


#####BEST MODEL Includes BOTH TUMI at lat 45 and some kind of climatic indicator

fit.2  <- gam(PHE.mean~MEI.YEAR.lag+TUMI.45.lag+Location.2, data=dat1)
summary(fit.2)
fit.2  <- gam(GLY.mean~MEI.YEAR.lag+TUMI.45.lag+Location.2, data=dat1)
summary(fit.2)
b <- getViz(fit.1)

#plot(sm(b,1)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ="darkgoldenrod3",lwd=1) + 
# l_ciLine(colour = "darkgoldenrod3", linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

plot(pterm(b,1)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ="darkgoldenrod3",lwd=1) + 
  l_ciLine(colour = "darkgoldenrod3", linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
plot(pterm(b,2)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ="darkgoldenrod3",lwd=1) + 
  l_ciLine(colour = "darkgoldenrod3", linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
plot(pterm(b, 3)) + l_ciBar(colour = c("skyblue3","darkgoldenrod3","skyblue3","darkgoldenrod3","skyblue3"), lwd=1) + 
  l_fitPoints(colour=c("skyblue3","darkgoldenrod3","skyblue3","darkgoldenrod3","skyblue3"), cex=3)+
  l_points(shape = 19, size = 1, alpha = 0.1)


plot(pterm(loc, 1)) + l_ciBar(colour = c("skyblue3","darkgoldenrod3","skyblue3","darkgoldenrod3","skyblue3"), lwd=1) + 
  l_fitPoints(colour=c("skyblue3","darkgoldenrod3","skyblue3","darkgoldenrod3","skyblue3"), cex=3)+
  l_points(shape = 19, size = 1, alpha = 0.1)













dat2 <- subset(dat1, Location.2=="SE")

fit.1 <- gam(PHE.mean~AK4.lag*AK2.lag*AK3.lag, data=dat1)
summary(fit.1)
AIC(fit.1) #432.9861*****STILL GOOD


#####BEST MODEL Includes BOTH TUMI at lat 45 and some kind of climatic indicator

fit.2 <- gam(PHE.mean~MEI.YEAR.lag+TUMI.45.lag+Location.2, data=dat1)
summary(fit.2)
b <- getViz(fit.1)

plot(sm(b,1)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ="darkgoldenrod3",lwd=1) + 
  l_ciLine(colour = "darkgoldenrod3", linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

plot(pterm(b,1)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ="darkgoldenrod3",lwd=1) + 
  l_ciLine(colour = "darkgoldenrod3", linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
plot(pterm(b,2)) +l_ciPoly(fill = "grey95")+ l_fitLine(colour ="darkgoldenrod3",lwd=1) + 
  l_ciLine(colour = "darkgoldenrod3", linetype = 2) + l_rug() +  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
plot(pterm(b, 3)) + l_ciBar(colour = c("skyblue3","darkgoldenrod3"), lwd=1) + 
  l_fitPoints(colour=c("skyblue4","darkgoldenrod4"), cex=3)+
  l_points(shape = 19, size = 1, alpha = 0.1)

loc <- gam(PHE.mean~Location.2, data=dat1)
summary(loc)
anova(loc)

plot(pterm(loc, 1)) + l_ciBar(colour = c("skyblue3","darkgoldenrod3"), lwd=1) + 
  l_fitPoints(colour=c("skyblue4","darkgoldenrod4"), cex=3)+
  l_points(shape = 19, size = 1, alpha = 0.1)

















