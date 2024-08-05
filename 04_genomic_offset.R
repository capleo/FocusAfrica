# Gradient forest
#PRELIMINARY
wd<-"~/focus/cowpea/"
setwd(wd)
rm(list=ls())
options(stringsAsFactors = F)

#load libraries
library(gradientForest)
library(rgdal)
library(raster)
library(sf)
library(RStoolbox)
library(ggspatial)
library(ggsflabel)
library(RColorBrewer)
library(maps)

library(tidyverse)
#library(ggplot2)
library(ggfortify)
library(patchwork)
library(ggpattern)

# load data + GF function
load("metadata.beforeGF.pruned.set.Rdata")
load("function.GF.cowpea.pruned.Rdata")
# load boarders
load(file="output/adm.crop.Rdata")

# create a folder to store results
GF.pruned <-"output/GF_pruned"
dir.create(GF.pruned)


#arrange variables by importance
by.imp.mem <- names(importance(gf.mem))

#Plot some summaries
# pdf("output/GF_pruned/GF_Overall importance.pdf")
plot(gf.mem, plot.type="Overall.Importance")
# dev.off()

##############################################################
# GF_pruned
##############################################################
#table(out$AEZ31)>10
crop_aez <- subset(aez, AEZ33 == 1|AEZ33 == 2)

# now crop and mask all model envs
clipped.modelEnv<-crop(modelEnv, extent(crop_aez))
clipped.modelEnv<-mask(clipped.modelEnv, crop_aez)
modelEnv<-clipped.modelEnv

clipped.modelEnv<-crop(modelFutureEnv_245, extent(crop_aez))
clipped.modelEnv<-mask(clipped.modelEnv, crop_aez)
modelFutureEnv_245<-clipped.modelEnv

clipped.modelEnv<-crop(modelFutureEnv_585, extent(crop_aez))
clipped.modelEnv<-mask(clipped.modelEnv, crop_aez)
modelFutureEnv_585<-clipped.modelEnv

##############################################################
# MAKE PREDICTIONS
##############################################################
# Predict allelic turnover actual environment
clim.land <- raster::extract(modelEnv, 1:ncell(modelEnv), df = TRUE) #clim.layer.crop
clim.land <- na.omit(clim.land)
pred.mem <- predict(gf.mem, clim.land[,-1])

# For projections with CMIP5, rcp85 and year 2050
clim.land.245 <- raster::extract(modelFutureEnv_245, 1:ncell(modelFutureEnv_245), df = TRUE) #clim.layer.crop
clim.land.245 <- na.omit(clim.land.245)
proj.mem.245 <- predict(gf.mem, clim.land.245[,-1]) 

# For projections with CMIP5, rcp85 and year 2070
clim.land.585 <- raster::extract(modelFutureEnv_585, 1:ncell(modelFutureEnv_585), df = TRUE) 
clim.land.585 <- na.omit(clim.land.585)
proj.mem.585 <- predict(gf.mem, clim.land.585[,-1]) 

#############################################################################################################
#Plot allelic turnover
#############################################################################################################
##These predictions then need to be converted to a color scale for mapping. 
##One way is to use principal components analysis (PCA) on the predictions and use 
##the first three axes to define red, green, blue color scales, respectively. 
##After the values are defined in color space, they can be stacked and mapped.

PCs <- prcomp(pred.mem, center=T, scale=F) #For pred.mem
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 +a2 -a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
mask<-modelEnv$bio3 #Precipitation of Coldest Quarter
mask[]<-as.numeric(mask[]>0)
rastR <- rastG <- rastB <- mask
rastR[clim.land$ID] <- r
rastG[clim.land$ID] <- g
rastB[clim.land$ID] <- b
rgb.rast <- stack(rastR, rastG, rastB)

# #For projection rcp85 year 2050
PCs.proj.245 <- prcomp(proj.mem.245, center=T, scale.=F) #For proj.mem
a1.proj.245 <- PCs.proj.245$x[, 1]
a2.proj.245 <- PCs.proj.245$x[, 2]
a3.proj.245 <- PCs.proj.245$x[, 3]
r.proj.245 <- a1.proj.245 + a2.proj.245
g.proj.245 <- -a2.proj.245
b.proj.245 <- a3.proj.245 +a2.proj.245 -a1.proj.245
r.proj.245 <- (r.proj.245 - min(r.proj.245))/(max(r.proj.245) - min(r.proj.245)) * 255
g.proj.245 <- (g.proj.245 - min(g.proj.245))/(max(g.proj.245) - min(g.proj.245)) * 255
b.proj.245 <- (b.proj.245 - min(b.proj.245))/(max(b.proj.245) - min(b.proj.245)) * 255
mask.proj.245<-modelFutureEnv_245$bio3 #Precipitation of Coldest Quarter
mask.proj.245[]<-as.numeric(mask.proj.245[]>0)
rastR.proj.245 <- rastG.proj.245 <- rastB.proj.245 <- mask.proj.245
rastR.proj.245[clim.land.245$ID] <- r.proj.245
rastG.proj.245[clim.land.245$ID] <- g.proj.245
rastB.proj.245[clim.land.245$ID] <- b.proj.245
rgb.rast.proj.245 <- stack(rastR.proj.245, rastG.proj.245, rastB.proj.245)

# #For projection rcp85 year 2070
PCs.proj.585 <- prcomp(proj.mem.585, center=T, scale.=F) #For proj.mem
r.proj.585 <- PCs.proj.585$x[, 1]
g.proj.585 <- PCs.proj.585$x[, 2]
b.proj.585 <- PCs.proj.585$x[, 3]
r.proj.585 <- (r.proj.585 - min(r.proj.585))/(max(r.proj.585) - min(r.proj.585)) * 255
g.proj.585 <- (g.proj.585 - min(g.proj.585))/(max(g.proj.585) - min(g.proj.585)) * 255
b.proj.585 <- (b.proj.585 - min(b.proj.585))/(max(b.proj.585) - min(b.proj.585)) * 255
mask.proj.585<-modelFutureEnv_585$bio3 #Precipitation of Coldest Quarter
mask.proj.585[]<-as.numeric(mask.proj.585[]>0)
rastR.proj.585 <- rastG.proj.585 <- rastB.proj.585 <- mask.proj.585
rastR.proj.585[clim.land.585$ID] <- r.proj.585
rastG.proj.585[clim.land.585$ID] <- g.proj.585
rastB.proj.585[clim.land.585$ID] <- b.proj.585
rgb.rast.proj.585 <- stack(rastR.proj.585, rastG.proj.585, rastB.proj.585)


#############################################################################################################
# Estimate genomic offset, so-called "genomic vulnerability"
#############################################################################################################

temp.245 <- vector("numeric", length = nrow(proj.mem.245))
for (i in 1:ncol(proj.mem.245)) {
  temp.245 <- temp.245 + (proj.mem.245[,i]-pred.mem[,i])^2
}


temp.585 <- vector("numeric", length = nrow(proj.mem.585))
for (i in 1:ncol(proj.mem.585)) {
  temp.585 <- temp.585 + (proj.mem.585[,i]-pred.mem[,i])^2
}

##############################################################
GenVuln.245 <- data.frame(sqrt(temp.245))
GenVuln.585 <- data.frame(sqrt(temp.585))

GenVuln <- cbind(clim.land[,c(1)], GenVuln.245, GenVuln.585)

colnames(GenVuln)[1] <- "cell_ID"
colnames(GenVuln)[2] <- "SSP245 (2041-2060)"
colnames(GenVuln)[3] <- "SSP 585 (2041-2060)"
summary(GenVuln)

# assign coordinates to each of the pixels
clim.land2 <- raster::extract(modelEnv, 1:ncell(modelEnv), df = TRUE)
clim.land2 <- cbind(coordinates(modelEnv), clim.land2)
clim.land2 <- na.omit(clim.land2)
genVuln <- cbind(clim.land2[,c(1,2)], GenVuln)

# Make GenVuln rasters
coordinates(genVuln) <- ~ x + y
gridded(genVuln) <- TRUE

genVuln.rast.245 <- raster(genVuln, "SSP245 (2041-2060)")
names(genVuln.rast.245)<-"Offset"
crs(genVuln.rast.245) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
#plot(genVuln.rast.245)

genVuln.rast.585 <- raster(genVuln, "SSP 585 (2041-2060)")
names(genVuln.rast.585)<-"Offset"

crs(genVuln.rast.585) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

plot(genVuln.rast.585)


