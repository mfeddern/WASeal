######################################       Trophic Regression      #####################################
data <- read.csv("Data/Compiled/SI Data.csv")

summary(lm(GLU.mean~ALA.mean, data=data))
summary(lm(GLU.mean~ASP.mean, data=data))
summary(lm(GLU.mean~VAL.mean, data=data))

summary(lm(ALA.mean~ASP.mean, data=data))
summary(lm(ALA.mean~VAL.mean, data=data))

summary(lm(ASP.mean~VAL.mean, data=data))

######################################        C3 verse C4 plant contributions        #####################################

vander <- read.csv("Data/Clean/VanderZanden.csv")
Function.MixingModel <- function(dataframe) {
  
  mixing.output <- ((dataframe[,"Obs"]-dataframe[,"C4EM"])/(dataframe[,"C3EM"]-dataframe[,"C4EM"]))*100
  
  return(mixing.output)
  
}




coastal <- data.frame(Obs=na.omit(subset(data, Location.2=="Coastal")$d13C.s), 
                      C4EM= rep(-9.5325, length(na.omit(subset(data, Location.2=="Coastal")$d13C.s))),
                      C3EM= rep(-19.44, length(na.omit(subset(data, Location.2=="Coastal")$d13C.s))))
salishsea <- data.frame(Obs=na.omit(subset(data, Location.2=="Inland")$d13C.s), 
                        C4EM= rep(-9.5325, length(na.omit(subset(data, Location.2=="Inland")$d13C.s))),
                        C3EM= rep(-19.44, length(na.omit(subset(data, Location.2=="Inland")$d13C.s))))
all <- data.frame(Obs=na.omit(subset(data, Location.2=="Coastal"|Location.2=="Inland")$d13C.s), 
                             C4EM= rep(-9.5325, length(na.omit(subset(data, Location.2=="Inland"|Location.2=="Coastal")$d13C.s))),
                             C3EM= rep(-19.44, length(na.omit(subset(data, Location.2=="Inland"|Location.2=="Coastal")$d13C.s))))

all.TEF <- data.frame(Obs=na.omit(subset(data, Location.2=="Coastal"|Location.2=="Inland")$d13C.s)-(3), 
                  C4EM= rep(-9.5325, length(na.omit(subset(data, Location.2=="Inland"|Location.2=="Coastal")$d13C.s))),
                  C3EM= rep(-19.44, length(na.omit(subset(data, Location.2=="Inland"|Location.2=="Coastal")$d13C.s))))

names2 <- c("Mean", "SD")
Percent.C3 <- data.frame(Coastal=c(mean(Function.MixingModel(coastal)), sd(Function.MixingModel(coastal))),
                          SalishSea=c(mean(Function.MixingModel(salishsea)), sd(Function.MixingModel(salishsea))),
                           All=c(mean(Function.MixingModel(all)), sd(Function.MixingModel(all))), None=c(100,100),
                           row.names = names2)


Function.MixingModel2 <- Function.MixingModel <- function(dataframe) {
  
  mixing.output <- ((dataframe[,"Obs"]-dataframe[,"C3EM"])/(dataframe[,"C4EM"]-dataframe[,"C3EM"]))*100
  
  return(mixing.output)
  
}

Percent.C4 <- data.frame(Coastal=c(mean(Function.MixingModel2(coastal)), sd(Function.MixingModel2(coastal))),
                         SalishSea=c(mean(Function.MixingModel2(salishsea)), sd(Function.MixingModel2(salishsea))),
                         All=c(mean(Function.MixingModel2(all)), sd(Function.MixingModel2(all))),
                         None=c(0,0),
                         row.names = names2)
Function.MixingModel2(all.TEF)


vander<-data.frame(vander)
str(vander)
mean.vander<-colMeans(na.omit(vander))
beta<-matrix(ncol =length(mean.vander)+1)

for(i in 1:length(mean.vander)){
  beta[,i] <- as.numeric(mean.vander[i]-mean.vander['Phe'])
}
colnames(beta)<- c(names(vander), "Pro")
beta<-data.frame(beta)

beta<-beta %>% select(Glu,
                       Ala,
                       Ile,
                       Leu,
                      Asp,
                      Val,
                      Pro)

beta[1,7]=-7.7

data2 <-data %>% select(PHE.mean,
                        GLU.mean,
                        ALA.mean,
                        VAL.mean,
                        ASP.mean,
                        PRO.mean,
                             d13C.s, 
                             Year,
                            
                             Location.2,
                        Sample.ID,
                        d13C,
                        d15N,
              
                        years)
data2 <- data2[complete.cases(data2), ]

all2 <- data.frame(Obs=na.omit(data2$d13C.s), 
                  C4EM= rep(-9.5325, length(na.omit(data2$d13C.s))),
                  C3EM= rep(-19.44, length(na.omit(data2$d13C.s))))

all3 <- data.frame(Obs=data2$d13C.s, 
                   C4EM= rep(-9.5325, length(data2$d13C.s)),
                   C3EM= rep(-19.44, length(data2$d13C.s)))

plot(Function.MixingModel2(all2),data2$PHE.mean, col=data2$Location.2, pch=16, xlab="Percent C4", ylab="d15NPHE")
legend(93, 20,legend=c("EBS", "WA Coast", "Salish Sea", "SC", "SE"),col=c("#ED90A4", "#C0AB52","#28BBD7", "#4FBF85",  "#C699E7"),lty=1,  cex=0.8)

plot(Function.MixingModel(all2),data2$PHE.mean, col=data2$Location.2, pch=16)
legend(90, 20,legend=c("EBS", "WA Coast", "Salish Sea", "SC", "SE"),col=c("#ED90A4", "#C0AB52","#28BBD7", "#4FBF85",  "#C699E7"),lty=1,  cex=0.8)



summary(lm(Function.MixingModel2(all2)~data2$PHE.mean))
summary(lmer(Function.MixingModel2(all2)~data2$PHE.mean+(1+PHE.mean|Location.2), data=data2))


palette(c("#ED90A4", "#C0AB52","#28BBD7", "#4FBF85",  "#C699E7"))
                       
########################    TEF Data                 ############################
data <- read.csv("Data/Compiled/SI Data.csv")
ger <- read.csv("Data/Compiled/GermainData2.csv")
ger.sealAA <- ger[1:7,7:19]
TEF.CH<- read.csv("Data/Compiled/TEFCh.csv")
JN<- read.csv("Data/Compiled/Nielsen_Beta.csv")
Beta.JN <- JN[1:2,]
TEF.JN <- JN[3:4,]

C3.weighted<- Beta.JN[1,2:8]*Percent.C3[1,1]/100
coastal.C4.weighted<-data.frame(beta*Percent.C4[1,1]/100)
salishsea.C4.weighted<-data.frame(beta*Percent.C4[2,1]/100)
all.C4.weighted<-data.frame(beta*Percent.C4[1,3]/100)

all.C4.weighted<-data.frame(beta*Percent.C4[1,3]/100)


coastal.beta <- cbind(Glu=C3.weighted['GLU']+coastal.C4.weighted['Glu'],
                      Val=C3.weighted['VAL']+coastal.C4.weighted['Val'],
                      Ala=C3.weighted['ALA']+coastal.C4.weighted['Ala'],
                      Asp=C3.weighted['ASP']+coastal.C4.weighted['Asp'],
                      Pro=C3.weighted['PRO']+coastal.C4.weighted['Pro'])

salishsea.beta <- cbind(Glu=C3.weighted['GLU']+salishsea.C4.weighted['Glu'],
                      Val=C3.weighted['VAL']+salishsea.C4.weighted['Val'],
                      Ala=C3.weighted['ALA']+salishsea.C4.weighted['Ala'],
                      Asp=C3.weighted['ASP']+salishsea.C4.weighted['Asp'],
                      Pro=C3.weighted['PRO']+salishsea.C4.weighted['Pro'])

all.beta <- cbind(Glu=C3.weighted['GLU']+all.C4.weighted['Glu'],
                        Val=C3.weighted['VAL']+all.C4.weighted['Val'],
                        Ala=C3.weighted['ALA']+all.C4.weighted['Ala'],
                   Asp=C3.weighted['ASP']+all.C4.weighted['Asp'],
                  Pro=C3.weighted['PRO']+all.C4.weighted['Pro'])



ind.beta<-matrix(nrow=length(Function.MixingModel2(all3)), ncol=5)
for(i in 1:length(Function.MixingModel2(all3))) {
  for(j in 1:5)
    ind.beta[i,j]<-((Function.MixingModel2(all3)/100)[i]*c(beta[1], beta[2], beta[5], beta[6],beta[7])[[j]])+
      ((1-(Function.MixingModel2(all3)/100)[i])*as.numeric(c(Beta.JN[1,2:3], Beta.JN[1,5:7]))[[j]])
}  
colnames(ind.beta)<- c("beta.glu", "beta.ala", "beta.asp", "beta.val","beta.pro")
ind<-cbind(ind.beta, data2,percC3=(1-(Function.MixingModel2(all3)/100)), LocC3=rep(NA, length(data2$Year)))
ind2 <- mutate(data, LocC3 = ifelse(Location.2 == "Coastal", Percent.C3[1,1], Percent.C3[1,2]))  
ind2 <- mutate(ind2, Loc.ala = ifelse(Location.2 == "Coastal", coastal.beta$ALA, salishsea.beta$ALA))  
ind2 <- mutate(ind2, Loc.glu = ifelse(Location.2 == "Coastal", coastal.beta$GLU, salishsea.beta$GLU))  
ind2 <- mutate(ind2, Loc.val = ifelse(Location.2 == "Coastal", coastal.beta$VAL, salishsea.beta$VAL))  
ind2 <- mutate(ind2, Loc.asp = ifelse(Location.2 == "Coastal", coastal.beta$ASP, salishsea.beta$ASP))  
ind2 <- mutate(ind2, Loc.pro = ifelse(Location.2 == "Coastal", coastal.beta$PRO, salishsea.beta$PRO))  


find.numeric <- sapply(ger.sealAA, is.numeric)
meanAA <- colMeans(ger.sealAA[, find.numeric])
hs <- data.frame(meanAA)
herring <- ger[8,7:19]
find.numeric <- sapply(herring, is.numeric)
herring <- colMeans(herring[, find.numeric])
herring<-data.frame(herring)
phe<- herring['Phe']
TEF <- data.frame(((hs-hs[10,])-herring+herring[10,]))
TEF<- data.frame(t(TEF))

########################  TEF EQ1                ############################

TP<-data.frame(rep(NA, length(AA.mean)))
Calculating.TP.1  <- function(dataframe, AA.mean, AA, AA2, x){
  for(i in 1:length(AA.mean)){
    TP[i,]<- (((AA.mean[i]-data$PHE.mean[i])-(Beta.JN[1,AA]*(x/100)+beta[AA2]*((100-x)/100)))/TEF.JN[1,AA])+1
  }
  TP 
}


########################   TEF EQ2                ############################
TP<-data.frame(rep(NA, length(AA.mean)))

Calculating.TP.2  <- function(dataframe, AA.mean, AA, AA2, x) {
  
  for(i in 1:length(AA.mean)){
    TP[i,]<- (((AA.mean[i]-data$PHE.mean[i])-(Beta.JN[1,AA]*(x/100)+beta[AA2]*((100-x)/100))-TEF[AA2])/
                TEF.JN[1,AA])+2
    
  }
  TP  
}



########################     Calculating TEF EQ3                ############################
TP<-data.frame(rep(NA, length(data$AA.mean)))

Calculating.TP.3  <- function(dataframe, AA.mean, AA, AA2, x) {
  
  for(i in 1:length(AA.mean)){
    TP[i,]<- (((AA.mean[i]-data$PHE.mean[i])-(Beta.JN[1,AA]*(x/100)+beta[AA2]*((100-x)/100))-TEF.CH[1,AA])/
                TEF.JN[1,AA])+3
    
  }
  TP  
}


########################     Calculating TEF EQ4               ############################
TP<-data.frame(rep(NA, length(data$AA.mean)))

Calculating.TP.4<- function(dataframe, AA.mean, AA, AA2, x) {
  
  for(i in 1:length(AA.mean)){
    TP[i,]<- (((AA.mean[i]-data$PHE.mean[i])-(Beta.JN[1,AA]*(x/100)+beta[AA2]*((100-x)/100))-TEF[AA2]-TEF.CH[1,AA])/
                TEF.JN[1,AA])+3
    
  }
  TP  
}



TP.4.beta<- data.frame(ALA=Calculating.TP.4(data, data$ALA.mean, 'ALA', 'Ala',x=Percent.C3[1,3]),
                       GLU=Calculating.TP.4(data, data$GLU.mean, 'GLU', 'Glu', x=Percent.C3[1,3]),
                       VAL=Calculating.TP.4(data, data$VAL.mean, 'VAL', 'Val', x=Percent.C3[1,3]),
                       ASP=Calculating.TP.4(data, data$ASP.mean, 'ASP', 'Asp', x=Percent.C3[1,3]),
                       PRO=Calculating.TP.4(data, data$PRO.mean, 'PRO', 'Pro', x=Percent.C3[1,3]))
colnames(TP.4.beta)<- c("TP.ALA", "TP.GLU", "TP.VAL", "TP.ASP", "TP.PRO")


########################     Calculating TEF EQ4 BETA IND              ############################


TP<-data.frame(rep(NA, length(ind$AA.mean)))

Calculating.TP.4IND<- function(dataframe, AA.mean, AA, AA2, beta) {
  
  for(i in 1:length(AA.mean)){
    TP[i,]<- (((AA.mean[i]-data$PHE.mean[i])-(beta[i])-TEF[AA2]-TEF.CH[1,AA])/
                TEF.JN[1,AA])+3
    
  }
  TP  
}

Calculating.TP.2IND<- function(dataframe, AA.mean, AA, AA2, beta) {
  
  for(i in 1:length(AA.mean)){
    TP[i,]<- (((AA.mean[i]-data$PHE.mean[i])-(beta[i])-TEF.CH[1,AA])/
                TEF.JN[1,AA])+2
    
  }
  TP  
}

TP.4IND<- data.frame(ALA=Calculating.TP.4IND(ind, ind$ALA.mean, 'ALA', 'Ala',ind$beta.ala),
                     GLU=Calculating.TP.4IND(ind, ind$GLU.mean, 'GLU', 'Glu', ind$beta.glu),
                     VAL=Calculating.TP.4IND(ind, ind$VAL.mean, 'VAL', 'Val', ind$beta.val),
                     ASP=Calculating.TP.4IND(ind, ind$ASP.mean, 'ASP', 'Asp', ind$beta.asp),
                     PRO=Calculating.TP.4IND(ind, ind$PRO.mean, 'PRO', 'Pro', ind$beta.pro),
                     
                     ALA=Calculating.TP.2IND(ind, ind$ALA.mean, 'ALA', 'Ala',ind$beta.ala),
                     GLU=Calculating.TP.2IND(ind, ind$GLU.mean, 'GLU', 'Glu', ind$beta.glu),
                     VAL=Calculating.TP.2IND(ind, ind$VAL.mean, 'VAL', 'Val', ind$beta.val),
                     ASP=Calculating.TP.2IND(ind, ind$ASP.mean, 'ASP', 'Asp', ind$beta.asp),
                     PRO=Calculating.TP.2IND(ind, ind$PRO.mean, 'PRO', 'Pro', ind$beta.pro))
colnames(TP.4IND)<- c("TP.ALA", "TP.GLU", "TP.VAL", "TP.ASP","TP.PRO",
                      "TP.ALA2", "TP.GLU2", "TP.VAL2", "TP.ASP2", "TP.PRO2")
                     
TP.LOC<- data.frame(ALA=Calculating.TP.4IND(ind2, ind2$ALA.mean, 'ALA', 'Ala',ind2$Loc.ala),
                     GLU=Calculating.TP.4IND(ind2, ind2$GLU.mean, 'GLU', 'Glu', ind2$Loc.glu),
                     VAL=Calculating.TP.4IND(ind2, ind2$VAL.mean, 'VAL', 'Val', ind2$Loc.val),
                     ASP=Calculating.TP.4IND(ind2, ind2$ASP.mean, 'ASP', 'Asp', ind2$Loc.asp),
                    PRO=Calculating.TP.4IND(ind2, ind2$PRO.mean, 'PRO', 'Pro', ind2$Loc.pro), 
                    
                     ALA=Calculating.TP.2IND(ind2, ind2$ALA.mean, 'ALA', 'Ala',ind2$Loc.ala),
                     GLU=Calculating.TP.2IND(ind2, ind2$GLU.mean, 'GLU', 'Glu', ind2$Loc.glu),
                     VAL=Calculating.TP.2IND(ind2, ind2$VAL.mean, 'VAL', 'Val', ind2$Loc.val),
                     ASP=Calculating.TP.2IND(ind2, ind2$ASP.mean, 'ASP', 'Asp', ind2$Loc.asp),
                    PRO=Calculating.TP.2IND(ind2, ind2$PRO.mean, 'PRO', 'Pro', ind2$Loc.pro))
colnames(TP.LOC)<- c(
                      "TP.ALA4.L", "TP.GLU4.L", "TP.VAL4.L", "TP.ASP4.L","TP.PRO4.L",
                      "TP.ALA2.L", "TP.GLU2.L", "TP.VAL2.L", "TP.ASP2.L","TP.PRO2.L")


########################   Full Dataframe                     ############################





TP <- data.frame(ALA=Calculating.TP.1(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,4]),
                  GLU=Calculating.TP.1(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,4]),
                  VAL=Calculating.TP.1(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,4]),
                  ASP=Calculating.TP.1(data, data$ASP.mean, 'ASP', 'Asp', Percent.C3[1,4]),
                 PRO=Calculating.TP.1(data, data$PRO.mean, 'PRO', 'Pro', Percent.C3[1,4]),
                  ALA=Calculating.TP.1(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,3]),
                  GLU=Calculating.TP.1(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,3]),
                  VAL=Calculating.TP.1(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,3]),
                  ASP=Calculating.TP.1(data, data$ASP.mean, 'ASP', 'Asp', Percent.C3[1,3]),
                 PRO=Calculating.TP.1(data, data$PRO.mean, 'PRO', 'Pro', Percent.C3[1,3]),
                  ALA=Calculating.TP.2(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,4]),
                  GLU=Calculating.TP.2(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,4]),
                  VAL=Calculating.TP.2(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,4]),
                  ASP=Calculating.TP.2(data, data$ASP.mean, 'ASP', 'Asp', Percent.C3[1,4]),
                 PRO=Calculating.TP.2(data, data$PRO.mean, 'PRO', 'Pro', Percent.C3[1,4]),
                  ALA=Calculating.TP.2(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,3]),
                  GLU=Calculating.TP.2(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,3]),
                  VAL=Calculating.TP.2(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,3]),
                  ASP=Calculating.TP.2(data, data$ASP.mean, 'ASP', 'Asp', Percent.C3[1,3]),
                 PRO=Calculating.TP.2(data, data$PRO.mean, 'PRO', 'Pro', Percent.C3[1,3]),
                  ALA=Calculating.TP.3(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,4]),
                  GLU=Calculating.TP.3(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,4]),
                  VAL=Calculating.TP.3(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,4]),
                  ASP=Calculating.TP.3(data, data$ASP.mean, 'ASP', 'Asp', Percent.C3[1,4]),
                 PRO=Calculating.TP.3(data, data$PRO.mean, 'PRO', 'Pro', Percent.C3[1,4]),
                  ALA=Calculating.TP.3(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,3]),
                  GLU=Calculating.TP.3(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,3]),
                  VAL=Calculating.TP.3(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,3]),
                  ASP=Calculating.TP.3(data, data$ASP.mean, 'ASP', 'Asp', Percent.C3[1,3]),
                 PRO=Calculating.TP.3(data, data$PRO.mean, 'PRO', 'Pro', Percent.C3[1,3]),
                  ALA=Calculating.TP.4(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,4]),
                  GLU=Calculating.TP.4(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,4]),
                  VAL=Calculating.TP.4(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,4]),
                  ASP=Calculating.TP.4(data, data$GLU.mean, 'ASP', 'Asp', Percent.C3[1,4]),
                 PRO=Calculating.TP.4(data, data$GLU.mean, 'PRO', 'Pro', Percent.C3[1,4]),
                  ALA=Calculating.TP.4(data, data$ALA.mean, 'ALA', 'Ala', Percent.C3[1,3]),
                  GLU=Calculating.TP.4(data, data$GLU.mean, 'GLU', 'Glu',  Percent.C3[1,3]),
                  VAL=Calculating.TP.4(data, data$VAL.mean, 'VAL', 'Val', Percent.C3[1,3]),
                  ASP=Calculating.TP.4(data, data$ASP.mean, 'ASP', 'Asp', Percent.C3[1,3]),
                PRO=Calculating.TP.4(data, data$PRO.mean, 'PRO', 'Pro', Percent.C3[1,3]
                                                            
                                       )
                )

colnames(TP)<- c("TP.ALA1", "TP.GLU1", "TP.VAL1", "TP.ASP1","TP.PRO1","TP.ALA1.beta", "TP.GLU1.beta", "TP.VAL1.beta", "TP.ASP1.beta","TP.PRO1.beta",
                   "TP.ALA2", "TP.GLU2", "TP.VAL2", "TP.ASP2","TP.PRO2","TP.ALA2.beta", "TP.GLU2.beta", "TP.VAL2.beta", "TP.ASP2.beta","TP.PRO2.beta",
                   "TP.ALA3", "TP.GLU3", "TP.VAL3", "TP.ASP3","TP.PRO3","TP.ALA3.beta", "TP.GLU3.beta", "TP.VAL3.beta", "TP.ASP3.beta","TP.PRO3.beta",
                   "TP.ALA4", "TP.GLU4", "TP.VAL4", "TP.ASP4","TP.PRO4","TP.ALA4.beta", "TP.GLU4.beta", "TP.VAL4.beta", "TP.ASP4.beta","TP.PRO4.beta")


########################    TP Plots variable TEF parameterization                      ############################
color<- rep(c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7'), each=4)
par(mfrow=c(1,1),mar=c(3,3,4,2))
cex=0.85
pdf(file="Results/Figures/TEF.pdf", width=5, height=5)
par(mfrow=c(2,2),mar=c(4,3,4,2))
plot(density(na.omit(TP$TP.ALA1)), xlim=c(1,7), ylim=c(0,2),xlab='Trophic Position',col='#7EBA68', cex.lab=0.75,lwd=2,lty=3, bty='n',cex.main=0.85,
     main=expression(paste("1. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe']), " - ", beta['Aq']),
                                      "TDF"['Average'])," + 1")))
polygon(density(na.omit(TP$TP.ALA1)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU1)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP1)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL1)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO1)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)
polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.2), border='grey')
par(xpd=TRUE)
legend(4.85, 2,legend=c("Alanine", "Glutamic Acid", "Valine", "Aspartic Acid", "Proline","Ecologically", "Realistic"),
       col=c('#7EBA68', '#CCA65A','#00C1B2', "#6FB1E7",'#D494E1',"grey",""),lty=c(3,1,2,4,5,1,0),  cex=0.7, bty='n')
par(xpd=FALSE)
sum(na.omit(TP[,1:5]) > 3.5 & na.omit(TP[,1:5]) < 5)/sum(na.omit(TP[,1:5]) > 0 & na.omit(TP[,1:5]) < 8)
text(4.25,1.85,labels="0.15",cex=cex)


plot(density(na.omit(TP$TP.ALA2)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68', lty=3, lwd=2, cex.lab=0.75,bty='n', xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("2. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ", "TDF"['HS'], " - ", beta['Aq']),
                                      "TDF"['Average'])," + 2")))
polygon(density(na.omit(TP$TP.ALA2)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU2)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP2)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL2)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO2)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
sum(na.omit(TP[,11:15]) > 3.5 & na.omit(TP[,11:15]) < 5)/sum(na.omit(TP[,11:15]) > 0 & na.omit(TP[,11:15]) < 8)
text(4.25,1.85,labels="0.22",cex=cex)


plot(density(na.omit(TP$TP.ALA3)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68', lwd=2, lty=3,cex.lab=0.75, bty='n', xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("3. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ", "TDF"['Phy'], " - ", beta['Aq']),
                                      "TDF"['Average'])," + 2")))
polygon(density(na.omit(TP$TP.ALA3)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU3)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP3)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL3)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO3)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
sum(na.omit(TP[,21:25]) > 3.5 & na.omit(TP[,21:25]) < 5)/sum(na.omit(TP[,21:25]) > 0 & na.omit(TP[,21:25]) < 8)
text(4.25,1.85,labels="0.66",cex=cex)


plot(density(na.omit(TP$TP.ALA4)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68',lty=3, lwd=2, cex.lab=0.75,bty='n', xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("4. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ", "TDF"['Phy'], " - ", "TDF"['HS']," - ", beta['Aq']),
                                      "TDF"['Average'])," + 3")))
polygon(density(na.omit(TP$TP.ALA4)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU4)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP4)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL4)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO4)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
sum(na.omit(TP[,31:35]) > 3.5 & na.omit(TP[,31:35]) < 5)/sum(na.omit(TP[,31:35]) > 0 & na.omit(TP[,31:35]) < 8)
text(4.25,1.85,labels="0.28",cex=cex)

dev.off()

########################    TP Plots variable TEF and Beta parameterization                    ############################
par(mfrow=c(1,1),mar=c(4,3,4,2))
pdf(file="Results/Figures/TEF.beta.pdf", width=5, height=5)
par(mfrow=c(2,2),mar=c(4,3,4,2))
plot(density(na.omit(TP$TP.ALA1.beta)), xlim=c(1,7), ylim=c(0,2),xlab='Trophic Position',col='#7EBA68',cex.lab=0.75, lwd=2,lty=3, bty='n',cex.main=0.85,
     main=expression(paste("1. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe']), " - ", beta['W']),
                                      "TDF"['Average'])," + 1")))

polygon(density(na.omit(TP$TP.ALA1.beta)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU1.beta)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP1.beta)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL1.beta)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO1.beta)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.2), border='grey')
par(xpd=TRUE)
legend(4.85, 2,legend=c("Alanine", "Glutamic Acid", "Valine", "Aspartic Acid", "Proline","Ecologically", "Realistic"),
       col=c('#7EBA68', '#CCA65A','#00C1B2', "#6FB1E7",'#D494E1',"grey",""),lty=c(3,1,2,4,5,1,0),  cex=0.7, bty='n')
par(xpd=FALSE)
sum(na.omit(TP[,6:10]) > 3.5 & na.omit(TP[,6:10]) < 5)/sum(na.omit(TP[,6:10]) > 0 & na.omit(TP[,6:10]) < 8)
text(4.25,1.85,labels="0.65",cex=cex)


plot(density(na.omit(TP$TP.ALA2.beta)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68', cex.lab=0.75,lty=3, lwd=2, bty='n', xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("2. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ", "TDF"['HS'], " - ", beta['W']),
                                      "TDF"['Average'])," + 2")))
polygon(density(na.omit(TP$TP.ALA2.beta)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU2.beta)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP2.beta)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL2.beta)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO2.beta)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')

sum(na.omit(TP[,16:20]) > 3.5 & na.omit(TP[,16:20]) < 5)/sum(na.omit(TP[,16:20]) > 0 & na.omit(TP[,16:20]) < 8)
text(4.25,1.85,labels="0.76",cex=cex)



plot(density(na.omit(TP$TP.ALA3.beta)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68', lwd=2, lty=3, bty='n', cex.lab=0.75,xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("3. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ", "TDF"['Phy'], " - ", beta['W']),
                                      "TDF"['Average'])," + 2")))
polygon(density(na.omit(TP$TP.ALA3.beta)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU3.beta)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP3.beta)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL3.beta)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO3.beta)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
sum(na.omit(TP[,26:30]) > 3.5 & na.omit(TP[,26:30]) < 5)/sum(na.omit(TP[,26:30]) > 0 & na.omit(TP[,26:30]) < 8)
text(4.25,1.85,labels="0.48",cex=cex)



plot(density(na.omit(TP$TP.ALA4.beta)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68',lty=3, lwd=2, bty='n',cex.lab=0.75, xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("4. ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ", "TDF"['Phy'], " - ", "TDF"['HS']," - ", beta['W']),
                                      "TDF"['Average'])," + 3")))
polygon(density(na.omit(TP$TP.ALA4.beta)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP$TP.GLU4.beta)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP$TP.ASP4.beta)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP$TP.VAL4.beta)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP$TP.PRO4.beta)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
sum(na.omit(TP[,36:40]) > 3.5 & na.omit(TP[,36:40]) < 5)/sum(na.omit(TP[,36:40]) > 0 & na.omit(TP[,36:40]) < 8)
text(4.25,1.85,labels="0.80",cex=cex)

dev.off()

########################     PLOT IND BETA                      ############################

pdf(file="Results/Figures/TEF.betaIND.pdf", width=5, height=5)
par(mfrow=c(2,2),mar=c(4,3,4,2))

plot(density(na.omit(TP.4IND$TP.ALA2)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68',lty=3, lwd=2, bty='n',cex.lab=0.75, xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("5.2 ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ",  "TDF"['HS']," - ", beta['Ind']),
                                       "TDF"['Average'])," + 2")))
polygon(density(na.omit(TP.4IND$TP.ALA2)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP.4IND$TP.GLU2)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP.4IND$TP.ASP2)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP.4IND$TP.VAL2)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP.4IND$TP.PRO2)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
par(xpd=TRUE)
legend(4.85, 2,legend=c("Alanine", "Glutamic Acid", "Valine", "Aspartic Acid", "Proline","Ecologically", "Realistic"),
       col=c('#7EBA68', '#CCA65A','#00C1B2', "#6FB1E7",'#D494E1',"grey",""),lty=c(3,1,2,4,5,1,0),  cex=0.7, bty='n')
par(xpd=FALSE)
sum(na.omit(TP.4IND[6:10]) > 3.5 & na.omit(TP.4IND[6:10]) < 5)/sum(na.omit(TP.4IND[6:10]) > 0 & na.omit(TP.4IND[6:10]) < 8)
text(4.25,1.85,labels="0.57",cex=cex)


plot(density(na.omit(TP.4IND$TP.ALA)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68',lty=3, lwd=2, cex.lab=0.75,bty='n', xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("5.2 ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ",  "TDF"['HS']," - ", beta['Ind']),
                                       "TDF"['Average'])," + 2")))
polygon(density(na.omit(TP.4IND$TP.ALA)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP.4IND$TP.GLU)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP.4IND$TP.ASP)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP.4IND$TP.VAL)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
polygon(density(na.omit(TP.4IND$TP.PRO)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

sum(na.omit(TP.4IND[1:5]) > 3.5 & na.omit(TP.4IND[1:5]) < 5)/sum(na.omit(TP.4IND[1:5]) > 0 & na.omit(TP.4IND[1:5]) < 8)
text(4.25,1.85,labels="0.70",cex=cex)


plot(density(na.omit(TP.LOC$TP.ALA2.L)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68',lty=3, lwd=2,cex.lab=0.75, bty='n', xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("5.2 ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ",  "TDF"['HS']," - ", beta['Loc']),
                                       "TDF"['Average'])," + 2")))
polygon(density(na.omit(TP.LOC$TP.ALA2.L)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP.LOC$TP.GLU2.L)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP.LOC$TP.ASP2.L)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP.LOC$TP.VAL2.L)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP.LOC$TP.PRO2.L)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)

polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
sum(na.omit(TP.LOC[6:10]) > 3.5 & na.omit(TP.LOC[6:10]) < 5)/sum(na.omit(TP.LOC[6:10]) > 0 & na.omit(TP.LOC[6:10]) < 8)
text(4.25,1.85,labels="0.40",cex=cex)


plot(density(na.omit(TP.LOC$TP.ALA4.L)), xlim=c(1,7), ylim=c(0,2),col='#7EBA68',lty=3, lwd=2,cex.lab=0.75, bty='n', xlab='Trophic Position',cex.main=0.85,
     main=expression(paste("5.4 ",frac(paste(paste(delta^15, "N"['Tr'])," - ",paste(delta^15, "N"['Phe'])," - ", "TDF"['Phy'], " - ", "TDF"['HS']," - ", beta['Loc']),
                                       "TDF"['Average'])," + 3")), cex.main=0.75)
polygon(density(na.omit(TP.LOC$TP.ALA4.L)),col=alpha('#7EBA68',0.25),border='#7EBA68', lwd=2, lty=3)
polygon(density(na.omit(TP.LOC$TP.GLU4.L)),col=alpha('#CCA65A',0.25),border='#CCA65A', lwd=2, lty=1)
polygon(density(na.omit(TP.LOC$TP.ASP4.L)),col=alpha('#00C1B2',0.25),border='#00C1B2',lwd=2,lty=2)
polygon(density(na.omit(TP.LOC$TP.VAL4.L)),col=alpha('#6FB1E7',0.25),border='#6FB1E7', lwd=2, lty=4)
polygon(density(na.omit(TP.LOC$TP.PRO4.L)),col=alpha('#D494E1',0.25),border='#D494E1', lwd=2, lty=6)
polygon(x=c(3.5,5,5,3.5), y=c(0,0,2,2), col=alpha('grey',0.15), border='grey')
sum(na.omit(TP.LOC[1:5]) > 3.5 & na.omit(TP.LOC[1:5]) < 5)/sum(na.omit(TP.LOC[1:5]) > 0 & na.omit(TP.LOC[1:5]) < 8)
text(4.25,1.85,labels="0.52",cex=cex)


dev.off()
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
########################     Creating Hierarchical Dataset subset selection       ############################

x<-0
b<-1
e<-1
Function.eq <- function(dataframe,dat, b, e, x){
                X1<- cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex", "Length")],dataframe[1+x], 
                      AA=rep("ALA",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
                colnames(X1)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex",  "Length", "TP", "AA", "eq", "beta")
                
               X2<- cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex", "Length")],dataframe[2+x], 
                      AA=rep("GLU",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
               colnames(X2)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex",  "Length", "TP", "AA", "eq", "beta")
               
               X3<- cbind(dat[,c("years","Location.2", "Sample.ID","PHE.mean", "d13C", "d13C.s", "d15N","Sex", "Length")],dataframe[3+x], 
                      AA=rep("VAL",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
               colnames(X3)<- c("years","Location.2", "Sample.ID","PHE.mean", "d13C", "d13C.s", "d15N","Sex",  "Length", "TP", "AA", "eq", "beta")
               
               X4<-  cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex", "Length")],dataframe[4+x], 
                      AA=rep("ASP",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
               colnames(X4)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex",  "Length", "TP", "AA", "eq", "beta")
               
               X5<-  cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex", "Length")],dataframe[5+x], 
                           AA=rep("PRO",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
               colnames(X5)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N","Sex",  "Length", "TP", "AA", "eq", "beta")
               
               output<- rbind(X1, X2, X3, X4, X5)
  return(output)
}



data.hier<- rbind(Function.eq(TP,data, 0, 1, 0),
Function.eq(TP,data,1, 1, 5),
Function.eq(TP,data,0, 2, 10),
Function.eq(TP,data,1, 2, 15),
Function.eq(TP,data,0, 3, 20),
Function.eq(TP,data,1, 3, 25),
Function.eq(TP,data,0, 4, 30),
Function.eq(TP,data,1, 4, 35),

Function.eq(TP.LOC,data,4, 4,0),
Function.eq(TP.LOC,data,4, 2,5))


Function.eq2 <- function(dataframe,dat, b, e, x){
  X1<- cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N")],dataframe[1+x], 
             AA=rep("ALA",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
  colnames(X1)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N", "TP", "AA", "eq", "beta")
  
  X2<- cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N")],dataframe[2+x], 
             AA=rep("GLU",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
  colnames(X2)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N", "TP", "AA", "eq", "beta")
  
  X3<- cbind(dat[,c("years","Location.2", "Sample.ID","PHE.mean", "d13C", "d13C.s", "d15N")],dataframe[3+x], 
             AA=rep("VAL",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
  colnames(X3)<- c("years","Location.2", "Sample.ID","PHE.mean", "d13C", "d13C.s", "d15N", "TP", "AA", "eq", "beta")
  
  X4<-  cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N")],dataframe[4+x], 
              AA=rep("ASP",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
  colnames(X4)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N", "TP", "AA", "eq", "beta")
  
  X5<-  cbind(dat[,c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N")],dataframe[5+x], 
              AA=rep("PRO",length(dataframe[1])), eq=rep(e, length(dataframe[1])), beta=rep(b, length(dataframe[1])))
  colnames(X4)<- c("years","Location.2", "Sample.ID", "PHE.mean","d13C", "d13C.s", "d15N", "TP", "AA", "eq", "beta")
  output<- rbind(X1, X2, X3, X4,X5)
  return(output)
}
data.hier2<- rbind(Function.eq2(TP.4IND,data2,3, 4, 0),
Function.eq2(TP.4IND,data2,3, 2, 5),
Function.eq2(TP.4IND,data2,4, 4,10),
Function.eq2(TP.4IND,data2,4, 2,15))





########################     Creating Hierarchical Dataset       ############################

lag <- 1
data4 <-cbind(Year=data.hier$years-lag, data.hier)
#total3 <- left_join(data4, data.d13C, by="Sample.ID")
#total3 <- left_join(total3,data.hPHE, by="Sample.ID")
#total3 <- left_join(total3,data.hd15N, by="Sample.ID")

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data4, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData.csv')




lag <- 2
data3 <-cbind(Year=data.hier$years-lag, data.hier)
#total3 <- left_join(data3, data.d13C, by="Sample.ID")
#total3 <- left_join(total3,data.hPHE, by="Sample.ID")
#total3 <- left_join(total3,data.hd15N, by="Sample.ID")

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data3, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData2.csv')


lag <- 3
data3 <-cbind(Year=data.hier$years--lag, data.hier)
prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data3, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData3.csv')


lag <- 4
data3 <-cbind(Year=data.hier$years-lag, data.hier)
#total3 <- left_join(data3, data.d13C, by="Sample.ID")
#total3 <- left_join(total3,data.hPHE, by="Sample.ID")
#total3 <- left_join(total3,data.hd15N, by="Sample.ID")

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data3, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData4.csv')


lag <- 5
data3 <-cbind(Year=data.hier$years+lag, data.hier)
#total3 <- left_join(data3, data.d13C, by="Sample.ID")
#total3 <- left_join(total3,data.hPHE, by="Sample.ID")
#total3 <- left_join(total3,data.hd15N, by="Sample.ID")

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data3, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData5.csv')

lag <- 6
data3 <-cbind(Year=data.hier$years+lag, data.hier)
#total3 <- left_join(data3, data.d13C, by="Sample.ID")
#total3 <- left_join(total3,data.hPHE, by="Sample.ID")
#total3 <- left_join(total3,data.hd15N, by="Sample.ID")

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data3, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData5.csv')

lag <- 7
data3 <-cbind(Year=data.hier$years+lag, data.hier)
#total3 <- left_join(data3, data.d13C, by="Sample.ID")
#total3 <- left_join(total3,data.hPHE, by="Sample.ID")
#total3 <- left_join(total3,data.hd15N, by="Sample.ID")

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data3, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/HierarchicalData5.csv')




########################     Creating Hierarchical Dataset       ############################

data6<-  cbind(data[,c("years","Location.2", "Sample.ID", "Sex", "Length", "d13C.s",
                       "d13C", "d15N", "PHE.mean", "GLU.mean", "VAL.mean", "PRO.mean",
                        "ASP.mean")],TP)
              


lag <- 2
data6 <-cbind(Year=data6$years+lag, data6)

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data6, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/DataFull.csv')


data7<-  cbind(data[,c("years","Location.2", "Sample.ID", "Sex", "Length", "d13C.s",
                       "d13C", "d15N", "PHE.mean", "GLU.mean", "VAL.mean", "PRO.mean",
                       "ASP.mean")],TP)
lag <- 0
data7 <-cbind(Year=data6$years+lag, data7)

prey <- read.csv("Data/Compiled/WA.Prey.tot.csv")
Env <- read.csv("Data/Compiled/Washington.Environmental.Standardized.csv")
total1 <- left_join(data7, Env, by="Year")
total2 <- left_join(total1, prey, by="Year")

write.csv(total2, 'Data/Compiled/DataFull0.csv')
