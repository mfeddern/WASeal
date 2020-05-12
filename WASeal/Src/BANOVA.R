rm(list = ls())
library(rjags)
library(runjags)
library(rstan)
library(BANOVA)

data <- read.csv("Data/Compiled/HierarchicalData.csv")
data2<-data.frame(cbind(data$Sample.ID, data$TP, data$AA))[1:805,]
colnames(data2)<- c("Sample.ID", "TP", "AA")
res2 <- BANOVA.run(TP~1,~AA, model_name = 'Normal', 
data =data2, id = 'Sample.ID', iter = 1000, chains = 2)


data2<- data2[complete.cases(data2), ]
app_1 <- BANOVA.Normal(TP~AA, data2,
                       + data2$Sample.ID, burnin = 5000, 
                       sample = 1000, thin = 20)
summary(aov(data2$TP~data2$AA))

BAnova(res2)
class(data)
data2$TP
