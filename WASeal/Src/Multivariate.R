data <- read.csv("Data/Compiled/HierarchicalData.csv")
#install.packages("brms")
library(dplyr)
library(brms)
#install.packages("mvtnorm") 

data("BTdata", package = "MCMCglmm")
head(BTdata)
fit1 <- brm(
  mvbind(tarsus, back) ~ sex + hatchdate + (1|p|fosternest) + (1|q|dam),
  data = BTdata, chains = 2, cores = 2
)
summary(fit1)

data5 <-data3 %>% select(allSmolt, 
                         HakeBiomass,
                         Herring.Biomass,
                         Chinook, 
                         HarborSeal,
                         TP.ALA,
                         TP.ASP,
                        # TP.PRO,
                        # TP.VAL,
                         TP, #TP for GLU
                         Year,
                         Sample.ID, 
                         Location.2)

Prey <- data5[complete.cases(data5), ]


fit1 <- brm(
  mvbind(TP, TP.ALA, TP.ASP) ~ Location.2+(1|q|Sample.ID),
  data = Prey, chains = 2, cores = 2, iter = 2500, control = list(max_treedepth = 10.5),
   prior=prior)

fit1 <- add_criterion(fit1, "loo")
summary(fit1)


get_prior(mvbind(TP, TP.ALA, TP.ASP) ~ Location.2+(1|q|Sample.ID),
          data = Prey)



prior <- set_prior(
  "normal(4,2)",
  class = "b",
  coef = "",
  group = "",
  resp = "",
  dpar = "",
  nlpar = "",
  lb =3,
  ub = 5,
  check = TRUE
)















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


fit1 <- brm(
  mvbind(TP.norm, TP.ALA.n, TP.ASP.n #, TP.VAL.n, TP.PRO.n
  ) ~ HarborSeal + Location.2+(1|q|Sample.ID),
  data = Prey, chains = 2, cores = 2, iter = 3500)

fit1 <- add_criterion(fit1, "loo")
summary(fit1)

