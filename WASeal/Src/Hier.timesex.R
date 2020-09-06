library(lme4)
library(AICcmodavg)
library("lme4")
library("ggplot2")
library("googleVis")
library("stargazer")
library("sjPlot")
library("lme4")
library('AICcmodavg')
########Fitting HIerarchical with sex and length
data <- read.csv("Data/Compiled/HierarchicalData2.csv")

dataSum <- subset(data, AA=="Glu"|AA=="PRO"|AA=="ALA"|AA=="VAL")
dataSum <- subset(dataSum,TP>0)

dataSum <- subset(dataSum,Sex=="F"|Sex=="M"&Year>1929)
dataSum <-dataSum %>% select(TP,
                               Year,
                               AA,
                               Sample.ID,
                               Location.2,
                              Sex)



dataSum <- dataSum[complete.cases(dataSum), ]
min(dataSum$Year)
max(dataSum$Year)
length(subset(dataSum, AA=='Glu')$Year)


AICc(lmer(TP~Year+Location.2+(1|AA), data=dataSum))
AICc(lmer(TP~Year+Location.2+Sex+(1|AA), data=dataSum))
AICc(lmer(TP~Sex+Location.2+(1|AA), data=dataSum))
AICc(lmer(TP~Location.2+(1|AA), data=dataSum))
AICc(lmer(TP~Sex*Location.2+(1|AA), data=dataSum))
AICc(lmer(TP~Year+I(Year^2)+Location.2+Sex+(1|AA), data=dataSum))
AICc(lmer(TP~Year+I(Year^2)+Location.2+(1|AA), data=dataSum))



Best <- lmer(TP~Year+Location.2+(1|AA), data=dataSum)
summary(Best)

color<- c('black','#CCA65A','#7EBA68','#00C1B2','#6FB1E7','#D494E1')
col<-ifelse(dataSum$Location.2=='Inland','#CCA65A','#6FB1E7')
AA.labs<-c("Alanine", "Glutamic Acid", "Proline", "Valine")
names(AA.labs)<-c("ALA", "Glu", "PRO", "VAL")
dataSum$MLPredictions <- fitted(Best)

gg <- ggplot(dataSum, aes(x = Year, y = TP, group = AA)) +
  labs(y="Trophic Position")+
  geom_point( alpha=0.5, size=3, color=col) +
  geom_smooth(aes(y=MLPredictions), se=FALSE,color = "black") +
  facet_wrap(~AA, labeller=labeller(AA=AA.labs))+
  theme_bw()+
  labs(colour="Location.2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(gg)

 year <- dataSum$Year
x.tilde=1
n.sims<- 1000
coef.hat<- as.matrix(coef(Best)$AA)[2,]
y.tilde.1928<- rnorm(1000, coef.hat%*%c(1, 1928, x.tilde),sigma.y.hat)
quantile(y.tilde.1928, c(0.25,0.5,0.75))
sd(y.tilde)

n.sims<- 1000
coef.hat<- as.matrix(coef(Best)$AA)[2,]
y.tilde.1977<- rnorm(1000, coef.hat%*%c(1, 1977, x.tilde),sigma.y.hat)
quantile(y.tilde.1977, c(0.25,0.5,0.75))
sd(y.tilde)

n.sims<- 1000
coef.hat<- as.matrix(coef(Best)$AA)[2,]
y.tilde.2014<- rnorm(1000, coef.hat%*%c(1, 2014, x.tilde),sigma.y.hat)
quantile(y.tilde.2014, c(0.25,0.5,0.75))
sd(y.tilde)


Coastal.pred <- c(mean(y.tilde.1928), mean(y.tilde.1977), mean(y.tilde.2014))


x.tilde=2
n.sims<- 1000
coef.hat<- as.matrix(coef(Best)$AA)[2,]
y.tilde.1928<- rnorm(1000, coef.hat%*%c(1, 1928, x.tilde),sigma.y.hat)
quantile(y.tilde.1928, c(0.25,0.5,0.75))
sd(y.tilde)

n.sims<- 1000
coef.hat<- as.matrix(coef(Best)$AA)[2,]
y.tilde.1977<- rnorm(1000, coef.hat%*%c(1, 1977, x.tilde),sigma.y.hat)
quantile(y.tilde.1977, c(0.25,0.5,0.75))
sd(y.tilde)

n.sims<- 1000
coef.hat<- as.matrix(coef(Best)$AA)[2,]
y.tilde.2014<- rnorm(1000, coef.hat%*%c(1, 2014, x.tilde),sigma.y.hat)
quantile(y.tilde.2014, c(0.25,0.5,0.75))
sd(y.tilde)
Inland.pred <- c(mean(y.tilde.1928), mean(y.tilde.1977), mean(y.tilde.2014))


data.table <-data.frame(Year=c(1928, 1977, 2014), Salish.Sea=Inland.pred, Coastal=Coastal.pred)
data.table<- digits
data.table2 <- data.table %>%
  gt()%>%
  fmt_number( # A column (numeric data)
    columns = vars(Salish.Sea, Coastal), # What column variable? BOD$Time
    decimals = 2 # With two decimal places
  ) %>% 
  tab_style(cell_text(size="medium"),  locations = cells_body(columns = vars(Year))) %>%
  tab_style(cell_text(size="medium"),  locations = cells_body(columns = vars(Coastal))) %>%
  tab_style(cell_text(size="medium"),  locations = cells_body(columns = vars(Salish.Sea))) %>%
  
  tab_header(
    title = md("Predicted Trophic Position"),
    subtitle = html("&sigma;<sub>y</sub><sup>2</sup> = 0.7"))%>%
cols_label(
  Year = md("**Year**"),
  Salish.Sea = md("**Salish Sea**"),
  Coastal = md("**Coastal**"))
data.table2
gtsave(data.table2, "Results/Tables/Pred.Table.GLU.png")







library(sm)
attach(mtcars)
dat<- data.frame(Year=c(rep(1928, n.sims), rep(1977, n.sims), 
                  rep(2014, n.sims)), TP=c(y.tilde.1928,
           y.tilde.1977, y.tilde.2014))
sm.density.compare(dat$TP, dat$Year)



