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
############### Time Series #################


dataTS <-data %>% select(
  TP,
  Location.2,
  AA,
  Year)%>%
  mutate(AA2 =recode(AA,'GLU' = "b. Glutamic Acid",'ALA'="a. Alanine",
                          'ASP'="d. Aspartic Acid", 'VAL'="c. Valine", "PRO"="e. Proline"))


dataTS <- dataTS[complete.cases(dataTS), ]
dataTS.ss<- subset(dataTS, Location.2=="Inland")
dataTS.ss$AA2<- factor(dataTS.ss$AA2, levels = c("b. Glutamic Acid","a. Alanine",
                                                 "d. Aspartic Acid", "c. Valine", "e. Proline"),
                              labels = c("Glutamic Acid","Alanine",
                                         "Aspartic Acid", "Valine", "Proline")
)
palette(c('#7EBA68', '#CCA65A','#00C1B2', "#6FB1E7",'#D494E1'))
palette()
col<-c('#7EBA68', '#CCA65A','#00C1B2', "#6FB1E7",'#D494E1')

dataTS.ss$Year <- as.numeric(dataTS.ss$Year)
SS.TS <- ggplot(dataTS.ss, aes(x =na.omit(Year), y = na.omit(TP))) + 
  stat_smooth(method = "gam", formula = TP~s(Year), data=dataTS.ss,aes(color = AA2)) +
  geom_point(aes(color = AA2, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA2~.,labeller = labeller(dataTS.ss$AA2))+
  ylim(1, 6)+
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
dataTS.c$AA2<- factor(dataTS.c$AA2,  levels = c("b. Glutamic Acid","a. Alanine",
                                                "d. Aspartic Acid", "c. Valine", "e. Proline"),
                      labels = c("Glutamic Acid","Alanine",
                                 "Aspartic Acid", "Valine", "Proline")
)
palette(c('#7EBA68', '#CCA65A','#00C1B2', "#6FB1E7",'#D494E1'))
palette()
Coastal.TS<-ggplot(dataTS.c, aes(x = Year, y = TP, color=col)) + 
  stat_smooth(fullrange=TRUE,aes(color = AA2)) +
  geom_point(aes(color = AA2, alpha=0.5), pch=16, size=2.5) +
  facet_grid(AA2~.,labeller = labeller(dataTS.c$AA2))+
  ylim(1, 6)+
  xlim(1925, 2015)+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = col)+
  ylab(" ")+
  ggtitle("B. Coastal")+
  guides(colour=FALSE, alpha=FALSE)
Coastal.TS

pdf(file="Results/Figures/TimeSeries.pdf", width=8, height=5.5)
ggarrange(SS.TS, Coastal.TS, 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 1, align="h")

dev.off()



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

summary(lm(TP~Sex, data=all.VAL.ss))
summary(lm(TP~Sex, data=all.GLU.ss))
summary(lm(TP~Sex, data=all.ALA.ss))
summary(lm(TP~Sex, data=all.ASP.ss))


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
                         eq)





Length.c <- Length%>% filter(Location.2 == "Coastal")
Length.c <- Length.c[complete.cases(Length.c), ]
Length.c$AA<- factor(Length.c$AA, levels = c("ALA","GLU",
                                             "ASP","VAL","PRO"),
                     labels = c("Alanine", "Glutamic Acid", "Aspartic Acid", "Valine","Proline")
)
summary(lm(TP~AA*Length, data=Length.c))
fit.c<-lm(TP~AA*Length, data=Length.c)

Coastal.Length <- qplot(Length, TP, data=Length.c, colour=AA)+
  ggtitle("A. Coastal")+
  theme_bw()+
  labs(y="Trophic Position", x="Standard Length (cm)")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  #geom_smooth(method="lm", aes(color = AA, alpha=0.5))+
  geom_line(data = fortify(fit.c), aes(x = Length, y = .fitted))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(Length.c$AA))+
  theme(plot.title = element_text(hjust = 0.5),strip.background =element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"
        ))+
  guides(colour=FALSE, alpha=FALSE)
Coastal.Length 


Length.ss <- Length%>% filter(Location.2 == "Inland")
Length.ss$AA<- factor(Length.ss$AA, levels = c("ALA","GLU",
                                               "ASP","VAL","PRO"),
                      labels = c("Alanine","Glutamic Acid",  "Aspartic Acid", "Valine", "Proline")
)
summary(lm(TP~AA*Length, data=Length.ss))
fit.ss<-lm(TP~AA*Length, data=Length.ss)

SalishSea.Length <- qplot(Length, TP, data=Length.ss, colour=AA)+
  ggtitle("B. Salish Sea")+
  theme_bw()+
  labs(y="", x="Standard Length (cm)")+
  geom_point(aes(color = AA, alpha=0.5), pch=16, size=2.5) +
  #geom_smooth(method="lm", aes(color = AA, alpha=0.5))+
  geom_line(data = fortify(fit.ss), aes(x = Length, y = .fitted))+
  scale_colour_manual(name = "AA",values = col)+
  facet_grid(AA~.,labeller = labeller(Length.ss$AA))+
  theme(plot.title = element_text(hjust = 0.5), strip.background =element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"
        ))+
  guides(colour=FALSE, alpha=FALSE)
SalishSea.Length 

pdf(file="Results/Figures/Lengthplot.pdf", width=8, height=4)
ggarrange(Coastal.Length, SalishSea.Length, rremove("x.text"), 
          #labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 1, align= 'hv')

dev.off()









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
