rm(list = ls()) 
data <- read.csv("Data/Processed/HerringSpawningBiomass.csv")
library(MARSS)
library(ggplot2)
library(forecast)

####Herring####
load(snotel.RData)
yy=data

acf(subset(data, Group=='CherryPoint'))
auto.arima(subset(data, Group=='CherryPoint')$Biomass)
acf(na.omit(subset(data, Group=='PortGamble')))
auto.arima(subset(data, Group=='PortGamble')$Biomass)
Arima(subset(data, Group=='PortGamble')$Biomass, order=c(1,0,0))

p <- ggplot(yy, aes(x = YEAR, y = Biomass)) + geom_line()
p + facet_wrap(~Group)

ns <- length(unique(yy$Group))
B <- "identity"
Q <- 'equalvarcov'
  #"equalvarcov"
R <- "diagonal and equal"
U <- "unequal"
A <- "zero"
x0 <- "unequal"

dat <- reshape2::acast(yy, Group ~ YEAR, value.var = "Biomass")


mod.list = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A)
dat.1 <- dat
dat.1[dat.1==0]<- NA

m <- apply(log(dat.1), 1, mean, na.rm = TRUE)
fit <- MARSS(log(dat.1), model = mod.list, control = list(maxit = 5000), 
             inits = list(A = matrix(m, ns, 1)))




library(broom)
library(ggplot2)
d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$Station <- d$.rownames
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                        ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")


Fit.1 <- fit
Fit.1$states.se
tot.herring.bio <- colSums(exp(Fit.1$states))
value<- tot.herring.bio
herring.tot <- cbind(yr=seq(1973,2012,1),value)
write.csv(herring.tot, "Data/Compiled/Prey/herring.tot.csv")
pdf(file="Results/Figures/HerringMarss.pdf", width=12, height=11)
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")
dev.off()

pdf(file="Results/Figures/HerringRaw.pdf", width=12, height=11)
p <- ggplot(yy, aes(x = YEAR, y = Biomass)) + geom_line()
p + facet_wrap(~Group)
dev.off()

####Harbor Seal####

rm(list = ls()) 
data <- read.csv("Data/Processed/HarborSeal.csv")
yy=data

pdf(file="Results/Figures/HarborSealRaw.pdf", width=12, height=11)
p <- ggplot(yy, aes(x = Year, y = Count)) + geom_line()
p + facet_wrap(~Group)
dev.off()

ns <- length(unique(yy$Group))
B <- "identity"
Q <- 'equalvarcov'
#"equalvarcov"
R <- "diagonal and equal"
U <- "unequal"
A <- "zero"
x0 <- "unequal"

dat <- reshape2::acast(yy, Group ~ Year, value.var = "Count")


mod.list = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A)
dat.1 <- dat
dat.1[dat.1==0]<- NA

m <- apply(log(dat.1), 1, mean, na.rm = TRUE)
fit <- MARSS(log(dat.1), model = mod.list, control = list(maxit = 5000), 
             inits = list(A = matrix(m, ns, 1)))




library(broom)
library(ggplot2)
d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$Station <- d$.rownames

pdf(file="Results/Figures/HarborSealMARSS.pdf", width=12, height=11)
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                          ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")
dev.off()


Fit.1 <- fit
Fit.1$states.se
tot.harborseal.pop <- colSums(exp(Fit.1$states))
value<- tot.harborseal.pop
seal.tot <- cbind(yr=seq(1975,1999,1),value)
write.csv(seal.tot, "Data/Compiled/Prey/seal.tot.csv")
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,5,1,1))
plot(seal.tot[,1], seal.tot[,2], ylab="Population", xlab="Year", type='l')


####Chinook####

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chinook2.csv")
yy <- cbind(Year=data[,1], Count =rowSums(data[,2:28]))
write.csv(yy, "Data/Compiled/Prey/chinook.total.csv")

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chinook3.csv")
yy=data

dat <- reshape2::acast(yy, Group ~ Year, value.var = "Spawners")


mod.list = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A)
dat.1 <- dat
dat.1[dat.1==0]<- NA

m <- apply(log(dat.1), 1, mean, na.rm = TRUE)
fit <- MARSS(log(dat.1), model = mod.list, control = list(maxit = 5000), 
             inits = list(A = matrix(m, ns, 1)))




####Chum####

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chum.csv")
yy <- cbind(Year=data[4:46,1], Count =rowSums(data[4:46,2:3]))
write.csv(yy, "Data/Compiled/Prey/chum.total.csv")

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chum2.csv")

pdf(file="Results/Figures/ChumRaw.pdf", width=12, height=11)
p <- ggplot(data, aes(x = Year, y = Spawners)) + geom_line()
p + facet_wrap(~Group)
dev.off()

####Coho####

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Coho.csv")


library(broom)
library(ggplot2)
d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$Station <- d$.rownames

pdf(file="Results/Figures/CohoMARSS.pdf", width=12, height=11)
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                          ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")
dev.off()

Fit.1 <- fit
Fit.1$states.se
tot.coho.pop <- colSums(exp(Fit.1$states))
value<- tot.coho.pop
coho.tot <- cbind(yr=seq(1957,2013,1),value)
write.csv(coho.tot, "Data/Compiled/Prey/coho.tot.csv")
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,5,1,1))
plot(coho.tot[,1], coho.tot[,2], ylab="Population", xlab="Year", type='l')

pdf(file="Results/Figures/CohoRaw.pdf", width=12, height=11)
p <- ggplot(data, aes(x = Year, y = Spawners)) + geom_line()
p + facet_wrap(~Group)
dev.off()

pdf(file="Results/Figures/CohoRaw.pdf", width=12, height=11)
p <- ggplot(data, aes(x = Year, y = Spawners)) + geom_line()
p + facet_wrap(~Group)
dev.off()
