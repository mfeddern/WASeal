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

#plotMap(nepacLL, xlim=c(-129, -115.6), ylim=c(43, 51.1),col="green",bg="lightblue")
color_palette<- function (n, h = c(360, -154), c = 100, l = c(43, 82), power = 1.5, 
                          fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
  if (!is.null(gamma)) 
    warning("'gamma' is deprecated and has no effect")
  if (n < 1L) 
    return(character(0L))
  h <- rep(h, length.out = 2L)
  c <- c[1L]
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, -1, length = n)
  rval <- hex(polarLUV(L = l[2L] - diff(l) * abs(rval)^power[2L], 
                       C = c * abs(rval)^power[1L], H = ifelse(rval > 0, h[1L], 
                                                               h[2L])), fixup = fixup, ...)
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}



col <- color_palette(7)

windows(width=5,height=7)  #you have to experiment with dimensions to get axes labels, etc.
#par(mfrow=c(1,1))

#plotMap(nepacLLhigh, xlim=c(-125, -121.9), ylim=c(46, 48.9),main="High res (nepacLLhigh)")
#compare worldHires in map package
#map("worldHires", xlim=c(-125, -121.9), ylim=c(46, 48.9),add=T,col="red")

#make sure to redraw map after resizing window!


#reading in shapefiles (uses maptools)
#water shapefile is from US Census (http://www.census.gov/geo/www/tiger/)
pdf(file="WAMapPoster.pdf", width=5, height=10)
mat <- matrix(data = c(1), ncol=1, nrow=1, byrow=TRUE)

lay <- layout(mat=mat)
#layout.show(3)
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
par(family = 'sans')

plotMap (nepacLLhigh, xlim=c(-125, -122), ylim=c(46, 49.5), col="white",bg=col[4], border='grey48', 
         axes = FALSE, ylab='', xlab='', colHoles='blue')


sanctuary <- importShapefile("ocnms_boundary/OCNMS_boundary.shp")

attributes(sanctuary)$projection <- "LL"

# prepare UTM coordinates matrix
utmcoor<-SpatialPoints(cbind(sanctuary$X,sanctuary$Y), proj4string=CRS("+proj=utm +zone=10"))
#utmdata$X and utmdata$Y are corresponding to UTM Easting and Northing, respectively.
#zone= UTM zone
# converting
longlatcoor<-data.frame(spTransform(utmcoor,CRS("+proj=longlat")))
sanc.latlong <- data.frame(cbind(sanctuary[,1],sanctuary[,2], sanctuary[,3], longlatcoor))
colnames(sanc.latlong) <- c("PID", "SID", "POS", "X", "Y")
attributes(sanc.latlong)$projection <- "LL"

#addPolys(sanc.latlong,col=col[4], border = "black", lty=1, lwd=1)

samples <- read.csv("Samples_InPosession_Summary.csv")

#samples.post <- subset(samples, Year >= 1978)
#samples.pre <- subset(samples, Year < 1978)

samples.i <- subset(samples, Location.2 == 'Inland')
samples.c <- subset(samples, Location.2 == 'Coastal')



t.1 <- adjustcolor(col[7], alpha.f=0.5)
t.2 <- adjustcolor('darkgoldenrod3', alpha.f=0.5)
t.14 <- adjustcolor(col[7], alpha.f=1)

#points(samples.pre$Long, samples.pre$Lat, pch=21, bg=paste(t.2), col=col[2])
#points(samples.post$Long, samples.post$Lat, pch=21, bg=paste(t.1), col=col[6])

points(samples.i$Long, samples.i$Lat, pch=16, col=paste(t.2), cex=2.25)
points(samples.c$Long, samples.c$Lat, pch=16, col=paste(t.1), cex=2.25)
#legend(-125.5,49.5, c('n = 1', expression("n ">="2"), "OCNMS"), 
  #      bg='white', col=c(paste(t.1), paste(t.14), 'black'), lty=c(0,0,1), pch=c(16,16,NA), pt.cex=1.5, cex=1.25)
dev.off()




coastal <- read.csv("Coastal.csv")
coastal[is.na(coastal)] <- 0

inland <- read.csv("Inland.csv")
inland[is.na(inland)] <- 0

par(oma=c(0,0,0,0), mar=c(5,5,4,1))
barplot(inland$Inland, ylab = "Number of Specimens", xlab="Year", col= 'darkgoldenrod', 
         ylim = c(0,12), main = "Salish Sea Specimens", cex.lab=1.5, 
        cex.axis = 1.5, cex.names = 1.5, cex.main=1.5)
axis(1, at=seq(0.7, 106.3, 1.2*4),labels=seq(1928, 2016, 4), pos=-0.5, cex.axis=1.25)


barplot(coastal$Coastal, ylab = "Number of Specimens", xlab="Year", col= col[7], 
         main = "Coastal Specimens", cex.lab=1.5, cex.axis = 1.5, 
        cex.names = 1.5, cex.main=1.5)
axis(1, at=seq(0.7, 106.3, 1.2*4),labels=seq(1928, 2016, 4), pos=-0.5, cex.axis=1.25)

#legend(175, 12, fill=c(col[2], col[5]),c('Salish Sea', 'Coastal'), col=c(col[5], col[2]






dev.off()

#library("colorspace") 
#pal <- choose_palette()
#pal(color_pa)











plot_ly(data, x = ~Year, y = ~Coastal, type = 'bar', name = 'Coastal', marker = list(color = col[2])) %>%
  add_trace(y = ~SalishSea, name = 'Salish Sea', marker = list(color = col[5])) %>%
  layout(yaxis = list(title = 'Count'), barmode = 'group')
