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
wildproduction <- read.csv("Data/Compiled/Prey/wildproduction.tot.csv")[2:3]

hatch2<- subset(hatch, Year>=1973&Year<=2013)
hatch2 <- cbind(hatch2$Smolts, wildproduction$value)
allsmolt <- cbind(wildproduction$yr, rowSums(hatch2))

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
wildproduction.norm <- cbind(Year=wildproduction[,1], Total=transform_to_log_scale(wildproduction[,2]-mean(na.omit(wildproduction[,2]))))
allsmolt.norm <- cbind(Year=allsmolt[,1], Total=transform_to_log_scale(allsmolt[,2]-mean(na.omit(allsmolt[,2]))))

herring.norm2 <- cbind(Year=herring$yr, Total=(herring$value-mean(herring$value))/sd(herring$value))
chinook.norm2 <- cbind(Year=chinook$Year, Total=(chinook$Count-mean(chinook$Count))/sd(chinook$Count))
chum.norm2 <- cbind(Year=chum$Year, Total=(chum$Count-mean(chum$Count))/sd(chum$Count))
coho.norm2 <- cbind(Year=coho$yr, Total=(coho$value-mean(coho$value))/sd(coho$value))
seal.norm2 <- cbind(Year=seal[,2], Total=(seal[,3]-mean(seal[,3]))/sd(seal[,3]))
hatch.norm2 <- cbind(Year=hatch$Year, Total=(hatch$Smolts-mean(hatch$Smolts))/sd(hatch$Smolts))
hake.norm2 <- cbind(Year=hake[,1], Total=(hake[,4]-mean(hake[,4]))/sd(hake[,4]))
eulachon.norm2 <- cbind(Year=eulachon[,1], Total=(eulachon[,2]-mean(eulachon[,2]))/sd(eulachon[,2]))
wildproduction.norm2 <- cbind(Year=wildproduction[,1], Total=(wildproduction[,2]-mean(wildproduction[,2]))/sd(wildproduction[,2]))
allsmolt.norm2 <- cbind(Year=allsmolt[,1], Total=(allsmolt[,2]-mean(allsmolt[,2]))/sd(allsmolt[,2]))



data.merge <- list(herring.norm, chinook.norm, chum.norm, coho.norm, seal.norm,
                   hatch.norm, hake.norm, eulachon.norm, wildproduction.norm, allsmolt.norm)
Wa.Prey <- Reduce(function(x,y) merge(x,y, by='Year', all.x=TRUE), data.merge)
colnames(Wa.Prey)<- c("Year","Herring.Biomass", "Chinook", "Chum", "Coho", "HarborSeal",
                      "HatcherySmolts", "HakeBiomass","EulachonLandings", "WildProduction", "allSmolt")

write.csv(Wa.Prey, "Data/Compiled/WA.Prey.tot.csv")


data.merge <- list(herring.norm2, chinook.norm2, chum.norm2, coho.norm2, seal.norm2,
                   hatch.norm2, hake.norm2, eulachon.norm2, wildproduction.norm2, allsmolt.norm2)
Wa.Prey2 <- Reduce(function(x,y) merge(x,y, by='Year', all.x=TRUE), data.merge)
colnames(Wa.Prey2)<- c("Year","Herring.Biomass", "Chinook", "Chum", "Coho", "HarborSeal",
                      "HatcherySmolts", "HakeBiomass","EulachonLandings", "WildProduction", "allSmolt")

write.csv(Wa.Prey2, "Data/Compiled/WA.Prey2.tot.csv")


