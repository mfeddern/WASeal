rm(list = ls())

#install.packages("mgcViz")
#install.packages("remotes")
#remotes::install_github("vapniks/mergeutils")
library(mgcViz)
library(AICcmodavg)
library(vapniks)
require(nlme)
require(mgcv)
library(qpcR)

############### Importing Datasets#########################
rm(list = ls())
lag <- 1

herring <- read.csv("Data/Compiled/Prey/herring.tot.csv")
chinook <- read.csv("Data/Compiled/Prey/chinook.total.csv")
chum <- read.csv("Data/Compiled/Prey/chum.total.csv")
seal <- read.csv("Data/Compiled/Prey/seal.tot.csv")
coho <- read.csv("Data/Compiled/Prey/coho.tot.csv")
hatch <- read.csv("Data/Compiled/Prey/salmon/RMIS_WA_smolts.csv")
hake <- read.csv("Data/Compiled/Prey/HakeNMFS.csv")
eulachon <- read.csv("Data/Compiled/Prey/Eulachon.csv")


########################## Normalize, demean############################

transform_to_log_scale <- function(x){
  if(x!=0){
    y <- (sign(x)) * (log(abs(x)))
  } else {
    y <-0
  }
  y 
}
herring.norm <- cbind(Year=herring$yr, Total=transform_to_log_scale(herring$value-mean(herring$value)))
chinook.norm <- cbind(Year=chinook$Year, Total=transform_to_log_scale(chinook$Count-mean(chinook$Count)))
chum.norm <- cbind(Year=chum$Year, Total=transform_to_log_scale(chum$Count-mean(chum$Count)))
coho.norm <- cbind(Year=coho$yr, Total=transform_to_log_scale(coho$value-mean(coho$value)))
seal.norm <- cbind(Year=seal[,2], Total=transform_to_log_scale(seal[,3]-mean(seal[,3])))
hatch.norm <- cbind(Year=hatch$Year, Total=transform_to_log_scale(hatch$Smolts-mean(na.omit(hatch$Smolts))))
hake.norm <- cbind(Year=hake[,1], Total=transform_to_log_scale(hake[,4]-mean(na.omit(hake[,4]))))
eulachon.norm <- cbind(Year=eulachon[,1], Total=transform_to_log_scale(eulachon[,2]-mean(na.omit(eulachon[,2]))))

data.merge <- list(herring.norm, chinook.norm, chum.norm, coho.norm, seal.norm, hatch.norm, hake.norm, eulachon.norm)
Wa.Prey <- Reduce(function(x,y) merge(x,y, by='Year', all.x=TRUE), data.merge)
colnames(Wa.Prey)<- c("Year","Herring.Biomass", "Chinook", "Chum", "Coho", "HarborSeal", "HatcherySmolts", "HakeBiomass","EulachonLandings")

write.csv(Wa.Prey, "Data/Compiled/WA.Prey.tot.csv")
