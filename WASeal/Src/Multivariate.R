data <- read.csv("Data/Compiled/HierarchicalData2.csv")
#install.packages("brms")
library(dplyr)
library(brms)
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
library(mvtnorm)
library(lme4)
library(tidyr)
library(dplyr)


data2<- subset(data, AA=='Glu')
#df <- tibble(x = 1:3, y = 3:1)
data3 <- cbind(data2, TP.ALA.n = subset(data, AA=='ALA')$TP.norm, 
               TP.VAL.n = subset(data, AA=='VAL')$TP.norm, 
               TP.ASP.n = subset(data, AA=='ASP')$TP.norm, 
               TP.PRO.n = subset(data, AA=='PRO')$TP.norm, 
               TP.ALA = subset(data, AA=='ALA')$TP, 
               TP.VAL = subset(data, AA=='VAL')$TP, 
               TP.ASP = subset(data, AA=='ASP')$TP, 
               TP.PRO = subset(data, AA=='PRO')$TP)
data5<- data3 %>% select(allSmolt, 
                         HakeBiomass,
                         Herring.Biomass,
                         Chinook, 
                         HarborSeal,
                         TP.ALA,
                        # TP.ASP,
                         TP.PRO,
                        #TP.VAL,
                         TP, #TP for GLU
                         Year,
                         Sample.ID, 
                         Location.2,
                       Chum,
                       Coho)

Prey <- data5[complete.cases(data5), ]



Prey.aic.output <- rbind(extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~allSmolt+Chinook+HakeBiomass+Location.2, data=Prey))[2],#7
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#8
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+allSmolt+Chinook+Location.2, data=Prey))[2],#9
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Herring.Biomass+Chinook+HakeBiomass+Location.2, data=Prey))[2],
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2+Chum, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2+Coho, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Coho+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Location.2+HarborSeal, data=Prey))[2], 
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Chinook+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+allSmolt+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+HakeBiomass+Location.2, data=Prey))[2],
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+Chinook+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+HakeBiomass+Location.2, data=Prey))[2],#2
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Herring.Biomass+allSmolt+Location.2, data=Prey))[2],#3
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+HakeBiomass+Chinook+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+allSmolt+Chinook+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+allSmolt+HakeBiomass+Location.2, data=Prey))[2],#6
                      
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Chum+Location.2, data=Prey))[2],#4
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~HarborSeal+Coho+Location.2, data=Prey))[2],#5
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Coho+Location.2, data=Prey))[2],#6
                      extractAIC(lm(cbind(TP, TP.ALA, TP.PRO) ~Chum+Coho+HarborSeal+Location.2, data=Prey))[2]#6
                      
                      
                      
  )
  
  model.names <- c("Location",
                   "Herring", "Chinook", "Hatch", "Hake","1. Herring, Chinook", "2. Herring, Hake", "3. Herring, Hatch", "4. Chinook, Hake",
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
  
  row.names(Prey.aic.output) <- model.names
  delaic <- Prey.aic.output-min(Prey.aic.output)
  aic.weight1 <- exp(-0.5*delaic)
  aic.weight <- aic.weight1/sum(aic.weight1)
  #dev.ex <- model.summary[,14]
  Prey.aic.output <- cbind(Prey.aic.output, delaic, aic.weight)
  
  colnames(Prey.aic.output)<- c("AIC", "delAIC", "AIC Weight")
  x<-data.frame(Prey.aic.output)
subset(x, delAIC<=2)

fit6 <- lm(cbind(TP, TP.ALA, TP.PRO) ~ Location.2+Herring.Biomass+HakeBiomass, data = Prey)
summary(fit6)
vcov(fit6)
Manova(fit6)
Anova(fit6)$aic


fit6 <- lm(cbind(TP, TP.ALA, TP.PRO) ~ Location.2, data = Prey)
summary(fit6)
vcov(fit6)
summary(Manova(fit6))
Anova(fit6)
extractAIC(fit6)


fit7 <- lm(cbind(TP, TP.ALA, TP.PRO) ~ Location.2+HakeBiomass, data = Prey)
extractAIC(fit7)


fit7 <- lmer(cbind(TP, TP.ALA, TP.PRO) ~ allSmolt+Location.2, data = Prey)
extractAIC(fit7)

fit1 <- add_criterion(fit1, "loo")
summary(fit6)

fit1 <- brm(mvbind(TP, TP.ALA, TP.PRO) ~ Location.2+HarborSeal+Herring.Biomass,
  data = Prey, chains = 1, cores = 2, iter = 2000)

get_prior(mvbind(TP, TP.ALA, TP.PRO) ~ Location.2+(1|q|Sample.ID),
          data = Prey)



#prior<- c(prior_string("normal(0,10)", class = "b", lb =3,ub = 5))

make_stancode(mvbind(TP, TP.ALA, TP.PRO) ~ Location.2+(1|q),
              data = Prey,
              prior = prior)











data2<- subset(data, AA=='Glu')
#df <- tibble(x = 1:3, y = 3:1)
data3 <- cbind(data2, TP.ALA.n = subset(data, AA=='ALA')$TP.norm, 
               TP.VAL.n = subset(data, AA=='VAL')$TP.norm, 
               TP.ASP.n = subset(data, AA=='ASP')$TP.norm, 
               TP.PRO.n = subset(data, AA=='PRO')$TP.norm, 
               TP.ALA = subset(data, AA=='ALA')$TP, 
               TP.VAL = subset(data, AA=='VAL')$TP, 
               TP.ASP = subset(data, AA=='ASP')$TP, 
               TP.PRO = subset(data, AA=='PRO')$TP)

data5 <-data3 %>% select(allSmolt, 
                         HakeBiomass,
                         Herring.Biomass,
                         Chinook, 
                         HarborSeal,
                         TP.ALA.n,
                         TP.ASP.n,
                         #TP.PRO.n,
                         #TP.VAL.n,
                         TP.norm,
                         Year,
                         Sample.ID,
                         Location.2)

Prey <- data5[complete.cases(data5), ]


fit2 <- brm(
  mvbind(TP, TP.ALA, TP.PRO #, TP.VAL.n, TP.PRO.n
  ) ~ Location.2+(1|q|Sample.ID),
  data = Prey, chains = 2, cores = 2, iter = 3500)

fit1 <- add_criterion(fit1, "loo")
summary(fit1)
loo(fit1,fit2)

fit1 <- brm(
  mvbind(TP, TP.ALA) ~ Location.2+(1|q|Sample.ID),
  data = data5, chains = 1, cores = 2, iter = 3500, 
  control = list(max_treedepth = 10.25, adapt_delta = 0.9))
fit2 <- brm(
  mvbind(TP, TP.ALA) ~ Location.2+Herring.Biomass+HarborSeal+(1|q|Sample.ID),
  data = Prey, chains = 1, cores = 2, iter = 3500, 
  control = list(max_treedepth = 10.25, adapt_delta = 0.9))
fit3 <- brm(
  mvbind(TP, TP.ALA) ~ Location.2+HarborSeal+(1|q|Sample.ID),
  data = Prey, chains = 1, cores = 2, iter = 3500, 
  control = list(max_treedepth = 10.25, adapt_delta = 0.9))
fit4 <- brm(
  mvbind(TP, TP.ALA) ~ Location.2+Herring.Biomass+(1|q|Sample.ID),
  data = Prey, chains = 1, cores = 2, iter = 3500, 
  control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo(fit1, fit2, fit3, fit4)
