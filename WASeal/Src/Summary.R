rm(list = ls())
install.packages("mgcv")
install.packages("ggplot2")
library(mgcViz)
library(AICcmodavg)
library(vapniks)
require(nlme)
require(mgcv)
library(qpcR)
library(dplyr)
library(dotwhisker)
library(ggpubr)

data <- read.csv("Data/Compiled/HierarchicalData.csv")
data<-subset(data, beta==1& eq==2)

############### Seasonality #################
data2 <- read.csv("Data/Compiled/WASealAAandTP2.csv")
data4 <- left_join(data, data2, by = "Sample.ID")

dataSeas <-data4 %>% select(
  TP,
  Location.2.x,
  AA,
  Year,
  Month)%>%
  mutate(AA2 =recode(AA,'GLU' = "a. Glutamic Acid",'ALA'="b. Alanine",
                     'ASP'="d. Aspartic Acid", 'VAL'="c. Valine", "PRO"="e. Proline"))

dataSeas$AA2<- factor(dataSeas$AA2, levels = c("a. Glutamic Acid","b. Alanine",
                                                 "d. Aspartic Acid", "c. Valine", "e. Proline"),
                       labels = c("Glutamic Acid","Alanine",
                                  "Aspartic Acid", "Valine", "Proline")
)

palette(c('#CCA65A','#7EBA68', '#00C1B2', "#6FB1E7",'#D494E1'))
palette()
col<-c('#CCA65A','#7EBA68', '#00C1B2', "#6FB1E7",'#D494E1')

dataSeas <- dataSeas[complete.cases(dataSeas), ]

Seas <- ggplot(dataSeas, aes(x=Month, y= TP, color=col)) + 
  #stat_smooth(aes(color = AA2), method='loess', formula=y~x)+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 12), aes(color = AA2))+
  geom_point(aes(color = AA2, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA2~.,labeller = labeller(dataSeas$AA2))+
  ylim(2, 7)+
  xlim(1, 12)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = col)+
  ylab("Trophic Position")+
  ggtitle("Seasonality")+
  guides(colour=FALSE, alpha=FALSE)
Seas

GLU<-filter(dataSeas, AA == "GLU")
summary(gam(TP~s(Month, k = 12), data=GLU))

ALA<-filter(dataSeas, AA == "ALA")
summary(gam(TP~s(Month, k = 12), data=ALA))

ASP<-filter(dataSeas, AA == "ASP")
summary(gam(TP~s(Month, k = 12), data=ASP))

VAL<-filter(dataSeas, AA == "VAL")
summary(gam(TP~s(Month, k = 12), data=VAL))

PRO<-filter(dataSeas, AA == "PRO")
summary(gam(TP~s(Month, k = 12), data=PRO))

pdf(file="Results/Figures/Seasonality.pdf", width=5, height=8)
Seas
dev.off()

############### Correlation Matrix #################
Corr.mat <-cbind(VAL=subset(data, AA=='VAL')$TP, PRO=subset(data, AA=='PRO')$TP,GLU=subset(data, AA=='GLU')$TP,
ALA=subset(data, AA=='ALA')$TP, ASP=subset(data, AA=='ASP')$TP)

cor(Corr.mat, method = "pearson", use = "complete.obs")
############### Time Series #################


dataTS <-data %>% select(
  TP,
  Location.2,
  AA,
  Year)%>%
  mutate(AA2 =recode(AA,'GLU' = "a. Glutamic Acid",'ALA'="b. Alanine",
                          'ASP'="d. Aspartic Acid", 'VAL'="c. Valine", "PRO"="e. Proline"))

dataTS<-subset(dataTS, Location.2=="Inland"|Location.2=="Coastal")
dataTS <- dataTS[complete.cases(dataTS), ]
dataTS.ss<- subset(dataTS, Location.2=="Inland")
dataTS.ss$AA2<- factor(dataTS.ss$AA2, levels = c("a. Glutamic Acid","b. Alanine",
                                                 "d. Aspartic Acid", "c. Valine", "e. Proline"),
                              labels = c("Glutamic Acid","Alanine",
                                         "Aspartic Acid", "Valine*", "Proline")
)
palette(c('#CCA65A','#7EBA68', '#00C1B2', "#6FB1E7",'#D494E1'))
palette()
col<-c('#CCA65A','#7EBA68', '#00C1B2', "#6FB1E7",'#D494E1')

dataTS.ss$Year <- as.numeric(dataTS.ss$Year)
gam.ss<- gam(TP~s(Year, k=5, by=AA), data=dataTS.ss)
ss.p<- predict_gam(gam.ss)

SS.TS <- ggplot(dataTS.ss, aes(x=Year, y= TP, color=col)) + 
  #stat_smooth(aes(color = AA2), method='loess', formula=y~x)+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 4), aes(color = AA2))+
  geom_point(aes(color = AA2, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA2~.,labeller = labeller(dataTS.ss$AA2))+
  ylim(2, 6)+
  xlim(1925, 2015)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = col)+
  ylab("Trophic Position")+
  ggtitle("A. Salish Sea")+
  guides(colour=FALSE, alpha=FALSE)
SS.TS

dataTS.c<- subset(dataTS, Location.2=="Coastal")
dataTS.c$AA2<- factor(dataTS.c$AA2,  levels = c("a. Glutamic Acid","b. Alanine",
                                                "d. Aspartic Acid", "c. Valine", "e. Proline"),
                      labels = c("Glutamic Acid*","Alanine",
                                 "Aspartic Acid*", "Valine", "Proline*")
)
palette(c('#7EBA68', '#CCA65A','#00C1B2', "#6FB1E7",'#D494E1'))
palette()
Coastal.TS<-ggplot(dataTS.c, aes(x = Year, y = TP, color=col)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 4), aes(color = AA2))+
  
  geom_point(aes(color = AA2, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA2~.,labeller = labeller(dataTS.c$AA2))+
  ylim(2, 6)+
  xlim(1925, 2015)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = col)+
  ylab(" ")+
  ggtitle("B. Coastal")+
  guides(colour=FALSE, alpha=FALSE)
Coastal.TS

pdf(file="Results/Figures/TimeSeries.pdf", width=8, height=8)
ggarrange(SS.TS, Coastal.TS, 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 1, align="h")

dev.off()

summary(gam(TP~s(Year, k=4), data=subset(dataTS.c, AA=="GLU")))#sig
summary(gam(TP~s(Year, k=4), data=subset(dataTS.c, AA=="ALA")))
summary(gam(TP~s(Year, k=4), data=subset(dataTS.c, AA=="ASP")))#sig
summary(gam(TP~s(Year, k=4), data=subset(dataTS.c, AA=="VAL")))
summary(gam(TP~s(Year, k=4), data=subset(dataTS.c, AA=="PRO")))#sig

summary(gam(TP~s(Year, k=4), data=subset(dataTS.ss, AA=="GLU")))
summary(gam(TP~s(Year, k=4), data=subset(dataTS.ss, AA=="ALA")))
summary(gam(TP~s(Year, k=4), data=subset(dataTS.ss, AA=="ASP")))
summary(gam(TP~s(Year, k=4), data=subset(dataTS.ss, AA=="VAL")))#sig
summary(gam(TP~s(Year, k=4), data=subset(dataTS.ss, AA=="PRO")))


Coastal.TS2<-ggplot(dataTS.c, aes(x = Year, y = TP, color=col)) + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 6), aes(color = AA2))+
  
  geom_point(aes(color = AA2, alpha=0.5), pch=16, size=2.5) +
  
  ylim(1, 6)+
  xlim(1925, 2015)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = col)+
  ylab(" ")+
  ggtitle("B. Coastal")+
  guides(colour=FALSE, alpha=FALSE)
Coastal.TS2
############### Location By Sex ###########
Sex<-data

males <- subset(subset(data, Sex=='M'), Location.2=="Inland"|Location.2=="Coastal")
length(males$Sex)/4 #54
females <- subset(subset(data, Sex=='F'), Location.2=="Inland"|Location.2=="Coastal")
length(females$Sex)/4 #61

col.f <-hcl.colors(6, palette = "Set 3", alpha = NULL, rev = FALSE, fixup = TRUE)
col.m <-hcl.colors(6, palette = "Dark 3", alpha = NULL, rev = FALSE, fixup = TRUE)

all.GLU.ss<- subset(Sex, AA=='GLU'& Location.2=="Inland")
male.GLU.ss <- subset(males, AA=='GLU'& Location.2=="Inland")
females.GLU.ss <- subset(females, AA=="GLU"& Location.2=="Inland")
all.ALA.ss<- subset(Sex, AA=='ALA'& Location.2=="Inland")
male.ALA.ss <- subset(males, AA=='ALA'& Location.2=="Inland")
females.ALA.ss <- subset(females, AA=="ALA"& Location.2=="Inland")
all.ASP.ss<- subset(Sex, AA=='ASP'& Location.2=="Inland")
male.ASP.ss <- subset(males, AA=='ASP'& Location.2=="Inland")
females.ASP.ss <- subset(females, AA=="ASP"& Location.2=="Inland")
all.VAL.ss<- subset(Sex, AA=='VAL'& Location.2=="Inland")
male.VAL.ss <- subset(males, AA=='VAL'& Location.2=="Inland")
females.VAL.ss <- subset(females, AA=="VAL"& Location.2=="Inland")
all.PRO.ss<- subset(Sex, AA=='PRO'& Location.2=="Inland")
male.PRO.ss <- subset(males, AA=='PRO'& Location.2=="Inland")
females.PRO.ss <- subset(females, AA=="PRO"& Location.2=="Inland")

all.GLU.c<- subset(Sex, AA=='GLU'& Location.2=="Coastal")
male.GLU.c <- subset(males, AA=='GLU'& Location.2=="Coastal")
females.GLU.c <- subset(females, AA=="GLU"& Location.2=="Coastal")
all.ALA.c<- subset(Sex, AA=='ALA'& Location.2=="Coastal")
male.ALA.c <- subset(males, AA=='ALA'& Location.2=="Coastal")
females.ALA.c <- subset(females, AA=="ALA"& Location.2=="Coastal")
all.ASP.c<- subset(Sex, AA=='ASP'& Location.2=="Coastal")
male.ASP.c <- subset(males, AA=='ASP'& Location.2=="Coastal")
females.ASP.c <- subset(females, AA=="ASP"& Location.2=="Coastal")
all.VAL.c<- subset(Sex, AA=='VAL'& Location.2=="Coastal")
male.VAL.c <- subset(males, AA=='VAL'& Location.2=="Coastal")
females.VAL.c <- subset(females, AA=="VAL"& Location.2=="Coastal")
all.PRO.c<- subset(Sex, AA=='PRO'& Location.2=="Coastal")
male.PRO.c <- subset(males, AA=='PRO'& Location.2=="Coastal")
females.PRO.c <- subset(females, AA=="PRO"& Location.2=="Coastal")


summary(lm(TP~Sex, data=all.VAL.c))
summary(lm(TP~Sex, data=all.GLU.c))
summary(lm(TP~Sex, data=all.ALA.c))
summary(lm(TP~Sex, data=all.ASP.c))
summary(lm(TP~Sex, data=all.PRO.c))

summary(lm(TP~Sex, data=all.VAL.ss))
summary(lm(TP~Sex, data=all.GLU.ss))
summary(lm(TP~Sex, data=all.ALA.ss))
summary(lm(TP~Sex, data=all.ASP.ss))
summary(lm(TP~Sex, data=all.PRO.ss))


n.sex <- length(male.GLU)+length(females.GLU) #102

Sex <-  subset(subset(subset(data, Sex=='F'|Sex=="M"), Location.2=="Inland"|Location.2=="Coastal"))
summary(lmer(TP~Sex+(1|AA), data=Sex))


pdf(file="Results/Figures/IsotopesbySexandLocation.pdf", width=10, height=8)

par(mfrow=c(2,1), mar=c(5,5,2,2))
boxplot(list(male.GLU.ss$TP, females.GLU.ss$TP,male.ALA.ss$TP, females.ALA.ss$TP,male.ASP.ss$TP, females.ASP.ss$TP,
             male.VAL.ss$TP, females.VAL.ss$TP, male.PRO.ss$TP, females.PRO.ss$TP),
        data=data, at = c(1,2, 4,5, 7,8, 10,11, 13, 14),
        main="A. Salish Sea",
        names = c("M", "F", "M", "F","M", "F","M", "F","M", "F"),
        col=(c(col.m[2], col.f[2], col.m[3], col.f[3], col.m[4], col.f[4],  col.m[5], col.f[5],  col.m[6], col.f[6])),
        pch=16, ylim=c(1, 6.5), ylab="Trophic Position")
text(1.5,1.5, labels = "Glutamic")
text(1.5,1, labels = "Acid")
text(4.5,1.5, labels = "Alanine")
#text(4.5,3, labels = "GoA")
text(7.5,1.5, labels = "Aspartic")
text(7.5,1, labels = "Acid")
text(10.5,1.5, labels = "Valine")
#text(1.5,20, labels = "*", cex=2)
#text(13.5,20, labels = "*", cex=2)
text(13.5,1.5, labels = "Proline")


boxplot(list(male.GLU.c$TP, females.GLU.c$TP,male.ALA.c$TP, females.ALA.c$TP,male.ASP.c$TP, females.ASP.c$TP,
             male.VAL.c$TP, females.VAL.c$TP,male.PRO.c$TP, females.PRO.c$TP),
        data=data, at = c(1,2, 4,5, 7,8, 10,11, 13, 14),
        names = c("M", "F", "M", "F","M", "F","M", "F","M", "F"),
        main="B. Coastal",
        col=(c(col.m[2], col.f[2], col.m[3], col.f[3], col.m[4], col.f[4],  col.m[5], col.f[5],  col.m[6], col.f[6])),
        pch=16, ylim=c(1, 6.5), ylab="Trophic Position")
text(1.5,1.5, labels = "Glutamic")
text(1.5,1, labels = "Acid")
text(4.5,1.5, labels = "Alanine")
#text(4.5,3, labels = "GoA")
text(7.5,1.5, labels = "Aspartic")
text(7.5,1, labels = "Acid")
text(10.5,1.5, labels = "Valine")
#text(1.5,20, labels = "*", cex=2)
#text(13.5,20, labels = "*", cex=2)
text(13.5,1.5, labels = "Proline")
dev.off()

Sex.ss<- subset(Sex, Location.2=="Inland")
summary(lmer(TP~Sex+(1|AA), data=Sex.ss))
AICc(lmer(TP~Sex+(1|AA), data=Sex.ss))
AICc(lmer(TP~(1|AA), data=Sex.ss))


Sex.c<- subset(Sex, Location.2=="Coastal")
summary(lmer(TP~Sex+(1|AA), data=Sex.c))
AICc(lmer(TP~Sex+(1|AA), data=Sex.ss))
AICc(lmer(TP~(1|AA), data=Sex.ss))
pairwise.t.test(data$TP, pool.sd=FALSE, data$AA, p.adj = "bonf")

############### Location By Length ###########
Length <-data %>% select(TP,
                         AA,
                         Length,
                         Location.2,
                         beta,
                         eq)%>%
                      mutate(AA2 =recode(AA,'GLU' = "a. Glutamic Acid",'ALA'="b. Alanine",
                     'ASP'="d. Aspartic Acid", 'VAL'="c. Valine", "PRO"="e. Proline"))

Length<- subset(Length, Location.2=="Inland"|Location.2=="Coastal")


Length.c <- Length%>% filter(Location.2 == "Coastal")
Length.c <- Length.c[complete.cases(Length.c), ]
Length.c$AA<- factor(Length.c$AA, levels = c("ALA","GLU",
                                             "ASP","VAL","PRO"),
                     labels = c("Alanine", "Glutamic Acid", "Aspartic Acid", "Valine","Proline")
)
summary(lm(TP~AA*Length, data=Length.c))
fit.c<-lm(TP~AA*Length, data=Length.c)

Coastal.Length <- ggplot(data=Length.c,aes(x=Length, y=TP, color=col))+
  ggtitle("B. Coastal")+
  theme_bw()+
  labs(y="Trophic Position", x="Standard Length (cm)")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  #geom_smooth(method="lm", aes(color = AA, alpha=0.5))+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 3), aes(color = AA))+
  #geom_line(data = fortify(fit.c), aes(x = Length, y = .fitted))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(Length.c$AA))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
Coastal.Length 


Length.ss <- Length%>% filter(Location.2 == "Inland")
Length.ss$AA<- factor(Length.ss$AA, levels = c("ALA","GLU",
                                               "ASP","VAL","PRO"),
                      labels = c("Alanine","Glutamic Acid",  "Aspartic Acid", "Valine", "Proline")
)
summary(lm(TP~AA*Length, data=Length.ss))
fit.ss<-lm(TP~AA*Length, data=Length.ss)

SalishSea.Length <- ggplot(data=Length.ss,aes(x=Length, y=TP, color=col))+
  ggtitle("A. Salish Sea")+
  theme_bw()+
  labs(y="Trophic Position", x="Standard Length (cm)")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  #geom_smooth(method="lm", aes(color = AA, alpha=0.5))+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 3), aes(color = AA))+
  #geom_line(data = fortify(fit.c), aes(x = Length, y = .fitted))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(Length.c$AA))+
  theme(plot.title = element_text(hjust = 0.5),strip.background =element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"
        ))+
  guides(colour=FALSE, alpha=FALSE)
SalishSea.Length 

Length.full <- Length

Length.full$AA2<- factor(Length.full$AA2, levels = c("a. Glutamic Acid","b. Alanine",
                                                 "d. Aspartic Acid", "c. Valine", "e. Proline"),
                       labels = c("Glutamic Acid","Alanine",
                                  "Aspartic Acid", "Valine*", "Proline")
)
summary(lm(TP~AA*Length, data=Length.ss))
fit.ss<-lm(TP~AA*Length, data=Length.ss)



Full.Length <- ggplot(data=Length.full,aes(x=Length, y=TP, color=col))+
  ggtitle("")+
  theme_bw()+
  labs(y="Trophic Position", x="Standard Length (cm)")+
  geom_point(aes(color = AA2, alpha=0.5), pch=16, size=2.5) +
  #geom_smooth(method="lm", aes(color = AA, alpha=0.5))+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 3), aes(color = AA2))+
  #geom_line(data = fortify(fit.c), aes(x = Length, y = .fitted))+
  scale_colour_manual(name = "AA2",values = col)+
  ylim(2,6)+
  facet_grid(AA2~.,labeller = labeller(Length.c$AA2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5))+
  guides(colour=FALSE, alpha=FALSE)
Full.Length 
pdf(file="Results/Figures/Lengthplot.pdf", width=4, height=6)
Full.Length
dev.off()


summary(gam(TP~s(Length, k=2), data=subset(Length.full, AA=="Glutamic Acid")))
summary(gam(TP~s(Length, k=4), data=subset(Length.full, AA=="Alanine")))
summary(gam(TP~s(Length, k=4), data=subset(Length.full, AA=="Aspartic Acid")))
summary(gam(TP~s(Length, k=4), data=subset(Length.full, AA=="Valine")))
summary(gam(TP~s(Length, k=4), data=subset(Length.full, AA=="Proline")))#sig

summary(lm(TP~log(Length), data=subset(Length.full, AA=="Glutamic Acid")))
summary(lm(TP~log(Length), data=subset(Length.full, AA=="Alanine")))
summary(lm(TP~log(Length), data=subset(Length.full, AA=="Aspartic Acid")))
summary(lm(TP~log(Length), data=subset(Length.full, AA=="Valine")))
summary(lm(TP~log(Length), data=subset(Length.full, AA=="Proline")))#sig




length.mod<-lmer(TP~Length+Location.2+(1|AA), data=subset(Length, Location.2=="Coastal"|Location.2=="Inland"))
Len
summary(lmer(TP~Length+Location.2+(1|AA), data=subset(Length, Location.2=="Coastal"|Location.2=="Inland")))
mean(subset(Length.full, Length >=150 & Length<=180 &AA=="Proline")$TP)#4.4
mean(subset(Length.full,  Length<=120 &AA=="Proline")$TP)#5.0
mean(subset(Length.c, Length >=150 & Length<=180 )$TP)#4.21
mean(subset(Length.c,  Length<=120 )$TP)#4.53


new.DATA.low <- data.frame(
  TP=Length$TP,
  Length = rep(120, length(Length$TP)), 
  Location.2=rep('Coastal', length(Length$TP)),  AA=Length$AA
)
pred.low <-data.frame(predict(length.mod, new.DATA.low, random.only=TRUE, interval = "confidence"))


new.DATA.high <- data.frame(
  TP=dataFull$TP,
  Col.Dis.high = dataFull$Col.Dis.high,
  allSmolt = rep(max(dataFull$allSmolt), 300), allSmolt=dataFull$allSmolt,
  Location.2=dataFull$Location.2, Col.Dis.high=dataFull$Col.Dis.high, AA=dataFull$AA
)
############### means and ss #################

mean(na.omit(subset(data, AA=="GLU")$TP))
sd(na.omit(subset(data, AA=="GLU")$TP))

mean(na.omit(subset(data, AA=="ALA")$TP))
sd(na.omit(subset(data, AA=="ALA")$TP))

mean(na.omit(subset(data, AA=="PRO")$TP))
sd(na.omit(subset(data, AA=="PRO")$TP))

mean(na.omit(subset(data, AA=="VAL")$TP))
sd(na.omit(subset(data, AA=="VAL")$TP))

mean(na.omit(subset(data, AA=="ASP")$TP))
sd(na.omit(subset(data, AA=="ASP")$TP))

############### Source AA verse TP AAs #################
data2 <- read.csv("Data/Compiled/DataFull.csv")
plot(data2$PHE.mean, data2$TP.GLU2.beta, ylim=c(0,6))
plot(data2$PHE.mean, data2$TP.ALA2.beta, ylim=c(0,6))
plot(data2$PHE.mean, data2$TP.PRO2.beta, ylim=c(0,6))
plot(data2$PHE.mean, data2$TP.ASP2.beta, ylim=c(0,6))
plot(data2$PHE.mean, data2$TP.VAL2.beta, ylim=c(0,6))


plot(data2$d15N, data2$TP.GLU2.beta, ylim=c(0,6), xlim=c(12,22))
plot(data2$d15N, data2$TP.ALA2.beta, ylim=c(0,6), xlim=c(12,22))
plot(data2$d15N, data2$TP.PRO2.beta, ylim=c(0,6), xlim=c(12,22))
plot(data2$d15N, data2$TP.ASP2.beta, ylim=c(0,6), xlim=c(12,22))
plot(data2$d15N, data2$TP.VAL2.beta, ylim=c(0,6), xlim=c(12,22))

plot(data2$PHE.mean, data2$GLU.mean)
plot(data2$PHE.mean, data2$ALA.mean)
plot(data2$PHE.mean, data2$PRO.mean, ylim=c(10,32))
plot(data2$PHE.mean, data2$VAL.mean)
plot(data2$PHE.mean, data2$ASP.mean, ylim=c(10,32))

plot(data2$GLU.mean, data2$TP.GLU2.beta, ylim=c(0,6),xlim=c(10,30))
plot(data2$ALA.mean, data2$TP.ALA2.beta, ylim=c(0,6),xlim=c(10,30))
plot(data2$PRO.mean, data2$TP.PRO2.beta, ylim=c(0,6),xlim=c(10,30))
plot(data2$ASP.mean, data2$TP.ASP2.beta, ylim=c(0,6),xlim=c(10,30))
plot(data2$VAL.mean, data2$TP.VAL2.beta, ylim=c(0,6),xlim=c(10,30))
