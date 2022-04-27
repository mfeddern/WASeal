### Packages ###
install.packages(c("maps","mapdata","maptools","mapproj"))
install.packages("PBSmapping")
install.packages("ggplot2")
install.packages("digest") #sometimes, I have to quit R, delete the digest directory from the library folder, then reinstall
devtools::install_github("dkahle/ggmap")
install.packages("RgoogleMaps")
install.packages('rgdal')
install.packages("colorspace")


library(maps)       #basic mapping functions and some data
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(rgdal)
library(colorspace)

#setwd("~/Documents/Graduate Work 2 copy/Classwork/Fall 2018/Map")
#setwd("~/Hardrive_Red/Graduate Work 2 copy/Classwork/Fall 2018/Map")
#######################################################################################
# PBS Mapping
#######################################################################################
install.packages("PBSmapping")
rm(list=ls())
graphics.off()

library(PBSmapping) #powerful mapping functions developed by Pacific Biological Station
data(nepacLL)
data(nepacLLhigh)


col <- c('#B88A00','#50A315','#00AD9A','#009ADE','#C86DD7')

pdf(file="Results/Figures/WAMap.pdf", width=14, height=9)
mat <- matrix(data = c(1,2,
                       1,3), ncol=2, nrow=2, byrow=TRUE)

lay <- layout(mat=mat)
#layout.show(3)
par(oma=c(0,0,0,0), mar=c(5,10,0,0))
par(family = 'sans')

plotMap(nepacLLhigh, xlim=c(-125.75, -121.5), ylim=c(46, 49.5),col="white", bg='grey90', cex.lab=1.5, cex.axis=1.35,
        border='grey48',  axes = FALSE, ylab='', xlab='', colHoles=NULL)
samples <- read.csv("Data/Compiled/Samples_InPosession_Summary.csv")

samples.i <- subset(samples, Location.2 == 'Inland')
samples.c <- subset(samples, Location.2 == 'Coastal')

t.1 <- adjustcolor(col[4], alpha.f=0.5)
t.2 <- adjustcolor(col[1], alpha.f=0.5)
t.14 <- adjustcolor(col[4], alpha.f=1)

#points(samples$Long, samples$Lat, pch=16, col=samples$Year, cex=2.5)


points(samples.i$Long, samples.i$Lat, pch=16, col=paste(t.2), cex=2.5)
points(samples.c$Long, samples.c$Lat, pch=16, col=paste(t.1), cex=2.5)
legend(-125.6,46.55, c('n = 1', expression("n ">="2"), "OCNMS"), 
        bg='white', col=c(paste(t.1), paste(t.14), 'black'), lty=c(0,0,1), pch=c(16,16,NA), pt.cex=1.5, cex=1.25)

coastal <- read.csv("Data/Compiled/MAPCoastal.csv")
coastal[is.na(coastal)] <- 0

inland <- read.csv("Data/Compiled/MAPInland.csv")
inland[is.na(inland)] <- 0

barplot(inland$Inland, ylab = "Number of Specimens", xlab="Year", col= col[1], 
        names.arg = coastal$Year, ylim = c(0,12), main = "Salish Sea Specimens", cex.lab=1.5, 
        cex.axis = 1.5, cex.names = 1.5, cex.main=1.5)

barplot(coastal$Coastal, ylab = "Number of Specimens", xlab="Year", col= col[4], 
        names.arg = coastal$Year, main = "Coastal Specimens", cex.lab=1.5, cex.axis = 1.5, 
        cex.names = 1.5, cex.main=1.5)
dev.off()




 #### SECOND MAP #####
samples <- read.csv("Data/Compiled/Samples_InPosession_Summary.csv")
samples <-samples %>% 
  mutate(Year.bin = cut(Year, breaks=c(1925, 1960, 1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2020)))
samples$Year.bin2<- factor(samples$Year.bin, levels = levels(samples$Year.bin),
                       labels = c("1925 - 1960", 
                                  "1960 - 1970", 
                                  "1970 - 1975", 
                                  "1975 - 1980", 
                                  "1980 - 1985",
                                  "1985 - 1990",
                                  "1990 - 1995",
                                  "1995 - 2000",
                                  "2000 - 2005",
                                  "2005 - 2020")
)
pdf(file="Results/Figures/WAMap2.pdf", width=4, height=4.5)
palette(hcl.colors(10,"Sunset", alpha=0.75))

par(oma=c(0,0,0,0), mar=c(5,10,0,0))
plotMap(nepacLLhigh, xlim=c(-127, -121.5), ylim=c(46, 49.5),col="white", bg='grey90', cex.lab=1.5, cex.axis=1.35,
        border='grey48',  axes = FALSE, ylab='', xlab='', colHoles=NULL)

points(samples$Long, samples$Lat, pch=16, col=as.factor(samples$Year.bin2), cex=1.5)

palette(hcl.colors(10,"Sunset", alpha=1))
legend(-127,48.55, pch=16, title= "Year Range", bg='white', legend=levels(samples$Year.bin2),col = as.factor(levels(samples$Year.bin2)))


dev.off()

