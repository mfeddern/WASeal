data <- read.csv("Data/Compiled/HierarchicalData.csv")
#install.packages("brms")
library(dplyr)
library(brms)
#install.packages("mvtnorm") 



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
                         TP.ALA,
                        # TP.ASP,
                         TP.PRO,
                        #TP.VAL,
                         TP, #TP for GLU
                         Year,
                         Sample.ID, 
                         Location.2)

#Prey <- data5[complete.cases(data5), ]

summary(lm(Prey$TP~Prey$TP.ALA))
summary(lm(Prey$TP~Prey$TP.ASP))
summary(lm(Prey$TP~Prey$TP.PRO))
summary(lm(Prey$TP~Prey$TP.VAL))
#Fit5, fit6, fit7, fit8

fit6 <- brm(
  mvbind(TP, TP.ALA, TP.PRO) ~ Location.2+Herring.Biomass+(1|q|Sample.ID),
  data = data5, chains = 1, cores = 2, iter = 3500, 
  control = list(max_treedepth = 10.25, adapt_delta = 0.9))



fit1 <- add_criterion(fit1, "loo")
summary(fit1)

fit1 <- brm(TP ~ Location.2+HarborSeal+Herring.Biomass+(1|q|Sample.ID),
  data = Prey, chains = 1, cores = 2, iter = 3500, 
  control = list(max_treedepth = 10.25, adapt_delta = 0.9))

get_prior(mvbind(TP, TP.ALA, TP.PRO) ~ Location.2+(1|q|Sample.ID),
          data = Prey)



#prior<- c(prior_string("normal(0,10)", class = "b", lb =3,ub = 5))

make_stancode(mvbind(TP, TP.ALA, TP.PRO) ~ Location.2+(1|q|Sample.ID),
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
