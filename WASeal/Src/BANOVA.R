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
data <- read.csv('Data/Compiled/Env.norm.csv')

period <- matrix(NA, nrow = length(data$TP.norm), ncol = 1)
data.p <-cbind(period, data)
subset(data.p, data.p$Year<=1945)

Period <- cbind(data.p, Period = dplyr::if_else(data.p['Year'] >=1995, data.p['period']<-1,
               if_else(data.p['Year'] <=1975 , data.p['period']<-3,
                       if_else(data.p['Year'] <=1948, data.p['period']<-4,
                       data.p['period']<-2))))



data4 <- Period.dat %>%
  select(-period) %>%
  select(-X) %>%
  select(-Sample.ID) %>%
  #select(-Location.2) %>%
  select(-Up.DFA1) %>%
  select(-Climate.DFA1.x) %>%
  select(-SST.DFA1)%>%
  select(-Dis.DFA1)%>%
  select(-years.x)%>%
  select(-years.y)%>%
  select(-PHE.norm)%>%
  select(-TP.norm)%>%
  select(-years)%>%
  select(-d13C.s)%>%
  select(-Year)%>%
  select(-PHE.mean)%>%
  select(-d13C.norm)

Wa.AOV <- data4[complete.cases(data4), ]
summary(aov(TP~ as.factor(Period) + Location.2+ Error(1/AA), data=Wa.AOV))

Wa.AOV2 <- subset(Wa.AOV, AA=='Glu'&Location.2=="Coastal")
summary(aov(TP~ as.factor(Period) + Location.2, data=Wa.AOV2))
summary(aov(TP~ as.factor(Period), data=Wa.AOV2))


mean(subset(Wa.AOV2, Period==1&Location.2=='Coastal')$TP)#4.44
mean(subset(Wa.AOV2, Period==2&Location.2=='Coastal')$TP)#5.72, consuming adult salmon? 
mean(subset(Wa.AOV2, Period==3&Location.2=='Coastal')$TP)#4.51
mean(subset(Wa.AOV2, Period==4&Location.2=='Coastal')$TP)#4.51


length(subset(Wa.AOV2, Period==2&Location.2=='Coastal')$TP)
length(subset(Wa.AOV2, Period==1&Location.2=='Coastal')$TP)
length(subset(Wa.AOV2, Period==3&Location.2=='Coastal')$TP)



mean(subset(Wa.AOV2, Period==1&Location.2=='Inland')$TP)#4.07
mean(subset(Wa.AOV2, Period==2&Location.2=='Inland')$TP)#4.00 
mean(subset(Wa.AOV2, Period==3&Location.2=='Inland')$TP)#3.79
mean(subset(Wa.AOV2, Period==4&Location.2=='Inland')$TP)#3.79

length(subset(Wa.AOV2, Period==2&Location.2=='Inland')$TP)
length(subset(Wa.AOV2, Period==1&Location.2=='Inland')$TP)
length(subset(Wa.AOV2, Period==3&Location.2=='Inland')$TP)




bf = anovaBF(TP~Location.2+AA, data=Wa.AOV, 
             whichRandom="AA")
Wa.AOV$AA
