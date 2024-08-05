# Gradient forest
#PRELIMINARY
wd<-"~/focus/cowpea/"
setwd(wd)
rm(list=ls())
options(stringsAsFactors = F)

#load libraries
library(raster)
library(maptools)
library(rgdal)

library(dismo)
library(data.table)

library(vegan)
library(BiodiversityR)
library(adespatial)
library(gradientForest)

library(ggplot2)
library(ggfortify)
library(patchwork)

#load data
load(file="output/metadata.pass.Rdata")

#load historical data
datafiles <- Sys.glob("input/bioclim/bioVars_Data/Historical/*.tif") #Or whatever identifies your files
currentEnv <- stack(datafiles)

datafiles_245 <- Sys.glob("input/bioclim/bioVars_Data/2041-2060/ssp245/*.tif") #Or whatever identifies your files
futureEnv_245 <- stack(datafiles_245)

datafiles_585<- Sys.glob("input/bioclim/bioVars_Data/2041-2060/ssp585/*.tif") #Or whatever identifies your files
futureEnv_585 <- stack(datafiles_585)

#get coordinates
coord <- coord[,c("LON","LAT")]

# estimate PCNMs
pcnm <- pcnm(dist(coord))  #this generates the PCNMs, you could stop here if you want all of them
keep <- ceiling(length(which(pcnm$value > 0))/2)
#keep <- round(length(which(pcnm$value > 0))/2)

pcnm.keep <- scores(pcnm)[,1:keep] #[,1:keep]    

# Calculate Moran's eigenvector maps (MEMs)
##This approach generates a set of uncorrelated spatial variables too

#mem <- dbmem(dist(coord))
mem2 <- dbmem(dist(coord), thresh=1.012, MEM.autocor = "positive", silent = FALSE)
names(attributes(mem2))

#barplot(attr(mem, "values"), main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)
barplot(attr(mem2, "values"), main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)

#Try to use same ratio of PCNMs... actually, it does not work
n.mem.keep <- ceiling(length(which(attr(mem2, "values") > 0))/2)

mem.keep <- scores(mem2)[,1:n.mem.keep] #keep half MEMs

###################################################################################################################
# now GRADIENT FOREST
###################################################################################################################
#we need same number of samples with bioclimatic variables and PCNMs or either MEMs
df.tmp <- meta[, c(2,12,22,23,3)]
#only georeferenced landraces
df.tmp <- na.omit(df.tmp)
# get rid of South African Acessions
df.tmp <- df.tmp[df.tmp$country != "South Africa", ] 

#rename the columns
colnames(df.tmp)[1] <- "ID"
colnames(df.tmp)[3] <- "LAT"
colnames(df.tmp)[4] <- "LON"

# exctract all historical climate data

pass<-SpatialPointsDataFrame(data=df.tmp, coords =df.tmp[,c("LON", "LAT")],
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  

# crop the rasters
modelEnv=crop(currentEnv,e)
modelFutureEnv_245=crop(futureEnv_245, e)
modelFutureEnv_585=crop(futureEnv_585, e)

##########################################
values <- raster::extract(modelEnv,pass)
passbio <- cbind(df.tmp[,c(1,3,4)],values)
##########################################

# Adjust names
names(modelEnv) <- sub("wc2p1_historical_","", names(modelEnv))
names(modelFutureEnv_245) <- sub("wc2p1_2041T2060_ssp245_","", names(modelFutureEnv_245))
names(modelFutureEnv_585) <- sub("wc2p1_2041T2060_ssp585_","", names(modelFutureEnv_585))

#create spatial points
sp<-SpatialPoints(coord)

bg <- randomPoints(modelEnv, nrow(coord))

##############################################################################################
#check redundancy in bioclim variables using the current values
modelEnv<-stack(modelEnv)

vif <- ensemble.VIF(
  x = modelEnv,
  a = data.frame(sp),
  VIF.max = 10,
  keep = NULL,
  layer.drops = NULL,
  factors = NULL,
  dummy.vars = NULL
  )

tokeep<-names(vif$VIF.final)

# drop colinear variables from environmental datasets
#redbc<-which(names(modelEnv) %in% tokeep)
modelEnv<-modelEnv[[tokeep]]
modelEnv<-brick(modelEnv)
modelFutureEnv_245<-modelFutureEnv_245[[tokeep]]
modelFutureEnv_585<-modelFutureEnv_585[[tokeep]]

colnames(passbio) <- gsub("wc2p1_historical_", "", colnames(passbio))

# create an object that contains only the climate and MEM spatial variables (no lat/lon)
env.gf.mem <- as.data.frame(cbind(passbio[,c("ID", tokeep)], mem.keep))

# check 
unique(env.gf.mem$ID[duplicated(env.gf.mem$ID)])

# first column as raw names
row.names(env.gf.mem) <- env.gf.mem[,1]
env.gf.mem <- env.gf.mem[,2:ncol(env.gf.mem)]

maxLevel.gf.mem <- log2(0.368*nrow(env.gf.mem)/2)

#Prepare a SNP file
snpGF <- read.table("input/snp.cowpea.forGF", header = T, row.names = 1)
snpGF[1:10,1:10]

#need 381 samples in the order of 
ordGF<-rownames(env.gf.mem)
snpGF<-subset(snpGF, rownames(snpGF) %in% ordGF)

env.gf.mem<-subset(env.gf.mem, rownames(env.gf.mem) %in% rownames(snpGF))

snpGF<-snpGF[order(match(rownames(snpGF), ordGF)), , drop = FALSE]
snpGF[1:20,1:10]
env.gf.mem[1:20,1:6]

# start with geno data, check missingness
perc.miss <- sum(is.na(snpGF))/(nrow(snpGF)*ncol(snpGF)) * 100

# simple, fast imputation (if total missiingness is minor than 5)
gen.imp <- apply(snpGF, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
gen.imp <- as.data.frame(gen.imp)

save.image(file="metadata.beforeGF.pruned.set.Rdata")

####################################################################################################
# RunGF
####################################################################################################
#server
setwd("/home/l.caproni/Documents/cowpea") #

load("metadata.beforeGF.pruned.set.Rdata")

gf.mem <- gradientForest(cbind(env.gf.mem, gen.imp),
                         predictor.vars=colnames(env.gf.mem),
                         response.vars=colnames(snpGF),
                         ntree=500,
                         maxLevel=maxLevel.gf.mem, trace=T,
                         corr.threshold=0.50,
                         nbin = 101, check.names = T)


save(gf.mem, file = "function.GF.cowpea.pruned.Rdata")

