rm(list = ls()) 
data <- read.csv("Data/qCov.csv")
library(MARSS)
library(ggplot2)
library(forecast)
head(data)
####Herring####
#load(snotel.RData)
fulldat =data

p <- ggplot(fulldat, aes(x = Time, y = Nitrate)) + geom_line()
p + facet_wrap(~River)

dat <-fulldat %>% select(Time,
                    Nitrate,
                    River)
dat <- reshape2::acast(dat, River ~ Time, value.var = "Nitrate")

covariates <-fulldat %>% select(River,
                                Qm3s,
                                Time)
covariates <- reshape2::acast(covariates, River ~ Time, value.var = "Qm3s")


the.mean = apply(dat, 1, mean, na.rm = TRUE)
the.sigma = sqrt(apply(dat, 1, var, na.rm = TRUE))
dat = (dat - the.mean) * (1/the.sigma)

the.mean = apply(covariates, 1, mean, na.rm = TRUE)
the.sigma = sqrt(apply(covariates, 1, var, na.rm = TRUE))
covariates = (covariates - the.mean) * (1/the.sigma)

d <- covariates
y <- dat  # to show relationship between dat & the equation
ns <- length(unique(yy$Group))
B <- "identity"
Q <- 'equalvarcov'
#"equalvarcov"
R <- "diagonal and equal"
U <- "unequal"
A <- "zero"
x0 <- "unequal"




model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, D = D, 
                   d = d, x0 = x0)
fit <- MARSS(y, model = model.list)

library(broom)
library(ggplot2)
d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$River <- d$.rownames
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                          ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~River) + xlab("") + ylab("SWE (demeaned)")




d <- covariates
y <- dat  # to show relationship between dat & the equation
ns <- length(unique(yy$Group))
B <- "identity"
Q <- 'equalvarcov'
#"equalvarcov"
R <- "diagonal and equal"
U <- "unequal"
A <- "zero"
x0 <- "unequal"


