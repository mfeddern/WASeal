data <- read.csv("Data/Compiled/HierarchicalData.csv")
#install.packages("brms")
library(dplyr)
library(brms)
library("bayesplot")
library("rstanarm")
library("ggplot2")
#install.packages("mvtnorm") 

########################    GLU PREY Models       ############################

dataPrey <- subset(data, Year>=1973&Year<=2008 & AA=="Glu")

data5 <-dataPrey %>% select(allSmolt, 
                         HakeBiomass,
                         Herring.Biomass,
                         Chinook, 
                         HarborSeal,
                         Coho,
                         Chum,
                         TP.norm,
                         TP,
                         Year,
                         Sample.ID,
                         Location.2)

Prey <- data5[complete.cases(data5), ]



data <- read.csv("Data/Compiled/HierarchicalData.csv")

AA <- brm(TP ~ AA-1, data = data, chains = 2, cores = 2, iter = 2000, 
         control = list(max_treedepth = 10.25, adapt_delta = 0.9))
posterior <- as.matrix(AA)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("b_AAALA", "b_AAASP", "b_AAGlu", "b_AAPRO", "b_AAVAL"),
           prob = 0.95) + plot_title

################# Candidate Models PREY ##########################################
L <- brm(TP ~ Location.2, data = Prey, chains = 2, cores = 2, iter = 2000, 
  control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.L<-loo(L)
He <- brm(TP ~ Location.2+Herring.Biomass,
            data = Prey, chains = 1, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.He<- loo(He)
C <- brm(TP ~ Chinook+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.C<- loo(C)


Sm <- brm(TP ~ allSmolt+Location.2,
           data = Prey, chains = 2, cores = 2, iter = 2000, 
           control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.Sm<- loo(Sm)
Ha <- brm(TP ~ HakeBiomass+Location.2,
           data = Prey, chains = 2, cores = 2, iter = 2000,
           control = list(max_treedepth = 10.25, adapt_delta = 0.9))##
loo.Ha<- loo(Ha)
He.C <- brm(TP ~ Herring.Biomass+Chinook+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.He.C<- loo(He.C)

He.Ha <- brm(TP ~Herring.Biomass+HakeBiomass+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))##
loo.He.Ha<- loo(He.Ha)

He.S <- brm(TP ~Herring.Biomass+allSmolt+Location.2,
             data = Prey, chains = 2, cores = 2, iter = 2000, 
             control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.He.S<- loo(He.S)

Ha.C <- brm(TP ~HakeBiomass+Chinook+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))##
loo.Ha.C<- loo(Ha.C)


S.C <- brm(TP ~allSmolt+Chinook+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.S.C<- loo(S.C)


S.Ha <- brm(TP ~allSmolt+HakeBiomass+Location.2,
           data = Prey, chains = 2, cores = 2, iter = 2000, 
           control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.S.Ha<- loo(S.Ha)


S.C.Ha <- brm(TP ~allSmolt+HakeBiomass+Chinook+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.S.C.Ha<- loo(S.C.Ha)


He.S.Ha <- brm(TP ~Herring.Biomass+allSmolt+HakeBiomass+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))##
loo.He.S.Ha<- loo(He.S.Ha)


He.S.C <- brm(TP ~Herring.Biomass+allSmolt+Chinook+Location.2,
               data = Prey, chains = 2, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.He.S.C<- loo(He.S.C)


He.Ha.C <- brm(TP ~Herring.Biomass+Chinook+HakeBiomass+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.He.Ha.C<- loo(He.Ha.C)


Hs.L <- brm(TP ~ HarborSeal+Location.2, data = data5, chains = 2, cores = 2, iter = 2000, 
         control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.Hs.L<- loo(Hs.L)


HS.He <- brm(TP ~ HarborSeal+Location.2+Herring.Biomass,
          data = Prey, chains = 2, cores = 2, iter = 2000, 
          control = list(max_treedepth = 10.25, adapt_delta = 0.9))##
loo.HS.He<- loo(HS.He)


HS.C <- brm(TP ~ HarborSeal+Chinook+Location.2,
         data = Prey, chains = 2, cores = 2, iter = 2000, 
         control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.C<- loo(HS.C)


HS.Sm <- brm(TP ~ HarborSeal+allSmolt+Location.2,
          data = Prey, chains = 2, cores = 2, iter = 2000, 
          control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.Sm<- loo(HS.Sm)


HS.Ha <- brm(TP ~ HarborSeal+HakeBiomass+Location.2,
          data = Prey, chains = 2, cores = 2, iter = 2000, 
          control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.Ha<- loo(HS.Ha)


HS.He.C <- brm(TP ~HarborSeal+ Herring.Biomass+Chinook+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.He.C<- loo(HS.He.C)


HS.He.Ha <- brm(TP ~HarborSeal+Herring.Biomass+HakeBiomass+Location.2,
             data = Prey, chains = 2, cores = 2, iter = 2000, 
             control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.He.Ha<- loo(HS.He.Ha)


HS.He.S <- brm(TP ~HarborSeal+Herring.Biomass+allSmolt+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.He.S<- loo(HS.He.S)


HS.Ha.C <- brm(TP ~HarborSeal+HakeBiomass+Chinook+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.Ha.C<- loo(HS.Ha.C)


HS.S.C <- brm(TP ~HarborSeal+allSmolt+Chinook+Location.2,
           data = Prey, chains = 2, cores = 2, iter = 2000, 
           control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.S.C<- loo(HS.S.C)


HS.S.Ha <- brm(TP ~HarborSeal+allSmolt+HakeBiomass+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.S.Ha<- loo(HS.S.Ha)


HS.S.C.Ha <- brm(TP ~HarborSeal+allSmolt+HakeBiomass+Chinook+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.S.C.Ha<- loo(HS.S.C.Ha)


HS.He.S.Ha <- brm(TP ~HarborSeal+Herring.Biomass+allSmolt+HakeBiomass+Location.2,
               data = Prey, chains = 2, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.He.S.Ha<- loo(HS.He.S.Ha)


HS.He.S.C <- brm(TP ~HarborSeal+Herring.Biomass+allSmolt+Chinook+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.He.S.C<- loo(HS.He.S.C)


HS.He.Ha.C <- brm(TP ~HarborSeal+Herring.Biomass+Chinook+HakeBiomass+Location.2,
               data = Prey, chains = 2, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo.HS.He.Ha.C<- loo(HS.He.Ha.C)


HS.He.Ha.C.S <- brm(TP ~HarborSeal+Herring.Biomass+Chinook+HakeBiomass+allSmolt+Location.2,
                  data = Prey, chains = 2, cores = 2, iter = 2000, 
                  control = list(max_treedepth = 10.25, adapt_delta = 0.9))

loo.HS.He.Ha.C.S<- loo(HS.He.Ha.C.S)

Ch <- brm(TP ~Chum+Location.2,
                    data = Prey, chains = 2, cores = 2, iter = 2000, 
                    control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.He <- brm(TP ~Chum+Herring.Biomass+Location.2,
          data = Prey, chains = 2, cores = 2, iter = 2000, 
          control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.C <- brm(TP ~Chum+Chinook+Location.2,
             data = Prey, chains = 2, cores = 2, iter = 2000, 
             control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.S <- brm(TP ~Chum+allSmolt+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.H <- brm(TP ~Chum+HakeBiomass+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.He.C <- brm(TP ~Chum+Herring.Biomass+Chinook+Location.2,
            data = Prey, chains = 2, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.He.H <- brm(TP ~Chum+Herring.Biomass+HakeBiomass+Location.2,
               data = Prey, chains = 2, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.He.S <- brm(TP ~Chum+allSmolt+Herring.Biomass+Location.2,
               data = Prey, chains = 2, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.H.C <- brm(TP ~Chum+allSmolt+HakeBiomass+Chinook+Location.2,
               data = Prey, chains = 2, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.S.C <- brm(TP ~Chum+allSmolt+Chinook+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))
Ch.S.H <- brm(TP ~Chum+allSmolt+HakeBiomass+Location.2,
              data = Prey, chains = 2, cores = 2, iter = 2000, 
              control = list(max_treedepth = 10.25, adapt_delta = 0.9))

HS.Ch <- brm(TP ~Chum+HarborSeal+Location.2,
          data = Prey, chains = 2, cores = 2, iter = 2000, 
          control = list(max_treedepth = 10.25, adapt_delta = 0.9))
HS.Ch.He <- brm(TP ~Chum+HarborSeal+Herring.Biomass+Location.2,
             data = Prey, chains = 2, cores = 2, iter = 2000, 
             control = list(max_treedepth = 10.25, adapt_delta = 0.9))
HS.Ch.C <- brm(TP ~Chum+HarborSeal+Chinook+Location.2,
            data = Prey, chains = 1, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
HS.Ch.S <- brm(TP ~Chum+HarborSeal+allSmolt+Location.2,
            data = Prey, chains = 1, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))
HS.Ch.H <- brm(TP ~Chum+HarborSeal+HakeBiomass+Location.2,
            data = Prey, chains = 1, cores = 2, iter = 2000, 
            control = list(max_treedepth = 10.25, adapt_delta = 0.9))

HS.He.C <- brm(TP ~Chum+HarborSeal+Herring.Biomass+Location.2,
                data = Prey, chains = 1, cores = 2, iter = 2000, 
                control = list(max_treedepth = 10.25, adapt_delta = 0.9))
HS.He.S <- brm(TP ~Chum+HarborSeal+Chinook+Location.2,
               data = Prey, chains = 1, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
HS.He.H <- brm(TP ~Chum+HarborSeal+allSmolt+Location.2,
               data = Prey, chains = 1, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))
HS.He <- brm(TP ~Chum+HarborSeal+HakeBiomass+Location.2,
               data = Prey, chains = 1, cores = 2, iter = 2000, 
               control = list(max_treedepth = 10.25, adapt_delta = 0.9))            


loo <- brms::loo(L, He, C, Sm, Ha, He.C, He.Ha, He.S, Ha.C, S.C, S.Ha, S.C.Ha, He.S.Ha, He.S.Ha,
    He.S.C, He.Ha.C, Hs.L, HS.He, HS.C, HS.Sm, HS.Ha, HS.He.C, HS.He.Ha, HS.He.S, 
    HS.Ha.C, HS.S.C, HS.S.Ha, HS.S.C.Ha, HS.He.S.Ha, HS.He.S.Ha, HS.He.S.C,HS.He.Ha.C,HS.He.Ha.C.S, reloo=TRUE)

loocomp <-loo_compare(loo.L, loo.He, loo.C, loo.Sm, loo.Ha, loo.He.C, loo.He.Ha, loo.He.S, loo.Ha.C, 
                      loo.S.C, loo.S.Ha, loo.S.C.Ha, loo.He.S.Ha, loo.He.S.Ha,
                      loo.He.S.C, loo.He.Ha.C, loo.Hs.L, loo.HS.He, loo.HS.C, loo.HS.Sm, loo.HS.Ha,
                      loo.HS.He.C, loo.HS.He.Ha, loo.HS.He.S, 
                      loo.HS.Ha.C,loo.HS.S.C, loo.HS.S.Ha, loo.HS.S.C.Ha, loo.HS.He.S.Ha, loo.HS.He.S.Ha,
                      loo.HS.He.S.C, loo.HS.He.Ha.C, loo.HS.He.Ha.C.S)



sig <- loocomp[,1]/loocomp[,2]
SE1<-subset(sig, sig > -1)
SE2<-subset(sig, sig > -2)




loo <- brms::loo(L, He, C,  Ha, He.C, He.Ha,  Ha.C, 
                 He.Ha.C, Hs.L, HS.He, HS.C, HS.Sm, HS.Ha, HS.He.C, HS.He.Ha, 
                 HS.Ha.C, HS.He.Ha.C, Ch, Ch.He, Ch.C,  Ch.H, Ch.He.C, Ch.He.H,
                 Ch.H.C, HS.Ch, HS.Ch.He, HS.Ch.C, HS.Ch.H,reloo=TRUE)



loo <- brms::loo(L, He, C, Sm, He.C, He.S, S.C,
                 He.S.C, Hs.L, HS.He, HS.C, HS.Sm,  HS.He.C,  HS.He.S, 
                  HS.S.C, HS.He.S.C, Ch, Ch.He, Ch.C, Ch.S, Ch.He.C, 
                 Ch.He.S,  HS.Ch, HS.Ch.He, HS.Ch.C, HS.Ch.S, reloo=TRUE)

# 
x <- loo_compare(loo.He, loo.C)
x[2,1]/x[2,2]
#Ch, Ch.He, Ch.C, Ch.S, Ch.H, Ch.He.C, Ch.He.H,
#Ch.He.S, Ch.H.C, Ch.S.H, HS.Ch, HS.Ch.He, HS.Ch.C, HS.Ch.S, HS.Ch.H,

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
                         TP,
                         TP.ALA,
                         TP.ASP,
                         TP.PRO,
                         TP.VAL,
                         TP.norm,
                         Year,
                         Sample.ID,
                         Location.2)

Prey <- data5[complete.cases(data5), ]


fit2 <- brm(
  mvbind(TP, TP.ALA, TP.PRO #, TP.VAL.n, TP.PRO.n
  ) ~ +Herring.Biomass+HarborSeal+(1|q|Sample.ID),
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


########################    GLU Nutrient Models       ############################

dataNutrient <- subset(data, Year>=1928&Year<=2014 & AA=="Glu")

data5 <-dataNutrient %>% select(d13C.norm,
                            d13C.s,
                            PHE.norm,
                            PHE.mean,
                            TP.norm,
                            TP,
                            Year,
                            Sample.ID,
                            Location.2)

Nutrient <- data5[complete.cases(data5), ]

################# Candidate Models Nutrient ##########################################
L2 <- brm(TP ~ Location.2, data = Nutrient, chains = 2, cores = 2, iter = 2000, 
         control = list(max_treedepth = 10.25, adapt_delta = 0.9))
P <- brm(TP ~ Location.2+PHE.norm, data = Nutrient, chains = 2, cores = 2, iter = 2000, 
         control = list(max_treedepth = 10.25, adapt_delta = 0.9))
C2 <- brm(TP ~ Location.2+d13C.norm, data = Nutrient, chains = 2, cores = 2, iter = 2000, 
         control = list(max_treedepth = 10.25, adapt_delta = 0.9))
CP <- brm(TP ~ Location.2+d13C.norm+PHE.norm, data = Nutrient, chains = 2, cores = 2, iter = 2000, 
          control = list(max_treedepth = 10.25, adapt_delta = 0.9))
loo(L2, P, C2, CP, reloo=TRUE)
