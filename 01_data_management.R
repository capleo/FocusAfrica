
# Prepare all geo data for subsequent anlysis

#PRELIMINARY
rm(list=ls())
wd<-"~/focus/cowpea/"
setwd(wd)

options(stringsAsFactors = F)

#load libraries
library(raster)
library(sf)
library(rgdal)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(ggrepel)
library(ggsflabel)

# load the shape of agroecologies
aez <- readOGR( 
  dsn= paste0(getwd(),"/input/AEZ33") , 
  layer="data_GAEZv4_33_GAEZv4_33_Africa.shp",
  verbose=FALSE
  )

#load labeles of Agroecologies

AEZ.lab <- read.delim (file="input/AEZ33_labeles.red.txt", header = T)

#load geo coordinates and passport data

meta <- read.csv(file="input/cowpea.info.pass_All.csv")

names(meta)
meta[1:10,1:10]

coord <- meta[, c(2,12,22,23,3)]

#only georeferenced landraces
coord <- na.omit(coord)

# get rid of South African Acessions
#table(coord$country)
coord <- coord[coord$country != "South Africa", ]      


#rename the columns
colnames(coord)[1] <- "ID"
colnames(coord)[3] <- "LAT"
colnames(coord)[4] <- "LON"


# crop the agroecologies to the study area the extent
e <- extent(min(coord$LON) -1, max(coord$LON) + 1,  min(coord$LAT) - 1 , max(coord$LAT) + 1)
aez.crop <- crop(aez, e)

# now convert the geo positions
pts<-st_as_sf(coord, coords = c("LON", "LAT"))
pts$LON<-coord$LON
pts$LAT<-coord$LAT

aez2 <- st_as_sf(aez.crop)
st_crs(pts)<-st_crs(aez2)

plot(pts)
head(pts)

#get intersection bw points and AEZs
inter<-st_intersection(aez2, pts)
head(inter)

#get relevant zones
dftmp<-data.frame(inter)
hitzones<-unique(dftmp[,"AEZ33"])

# add a factor for plotting
aez2 <- merge(aez2, AEZ.lab, by.x="AEZ33", by.y="AEZ")

# aez2$AEZ<-as.character(aez2$AEZ)
# aez2$AEZ[which(!aez2$AEZ %in% hitzones)]<-"ZNR"
aez2$AEZ<-as.character(aez2$AEZdef)
aez2$AEZ[which(!aez2$AEZ %in% AEZ.lab$AEZdef)]<-"ZNR"
aez2$AEZ<-as.factor(aez2$AEZ)

# colour palette for the AEZs
nicecols3<- c ("#d9d778" ,"#d0dfb2","#74e451","#d3dd3e", 
               "#657452", "#47c39d" , "#de9342")

# #make a nice map
# outmap<-ggplot(data = aez2) +
#   geom_sf(aes(fill=AEZ), alpha=0.8) + 
#   geom_sf(data=pts, size=1.8, color="red", alpha=0.6) +
#   labs(fill = "Agroecological \nZones", size = "") +
#   xlab("longitude") + ylab("latitude")+
#   theme_light() +
#   scale_fill_manual(values=c(nicecols3, "lightblue")) +
#   #theme(legend.position="bottom") +
#   #geom_sf(data=pts_pheno, shape=23, fill="red", color="black", size=2.5, alpha=0.7 ) +
#   #geom_sf_label_repel(data = pts_pheno, aes(label = loc),nudge_x = +5.0, nudge_y = 0, seed = 10, alpha =0.7)+
#   guides(color=guide_legend(title="DAPC cluster")) +
#   annotation_scale(location = "br", width_hint = 0.2) +
#   annotation_north_arrow(location = "br", which_north = "true", 
#                          pad_x = unit(0, "in"), pad_y = unit(0.2, "in"),
#                          style = north_arrow_fancy_orienteering) 
# outmap

pts$DAPC <- as.factor(pts$DAPC)

#make a nice map with clusters
outmap<-ggplot(data = aez2) +
  geom_sf(aes(fill=AEZ), alpha=0.8) + 
  geom_sf(data=pts, size=1.8, aes(color=DAPC), alpha=0.6) +
  labs(fill = "Agroecological \nZones", size = "") +
  xlab("longitude") + ylab("latitude")+
  theme_light() +
  scale_fill_manual(values=c(nicecols3, "lightblue")) +
  scale_color_manual(values = c("red", "green", "blue")) +
  #theme(legend.position="bottom") +
  #geom_sf(data=pts_pheno, shape=23, fill="red", color="black", size=2.5, alpha=0.7 ) +
  #geom_sf_label_repel(data = pts_pheno, aes(label = loc),nudge_x = +5.0, nudge_y = 0, seed = 10, alpha =0.7)+
  guides(color=guide_legend(title="DAPC cluster")) +
  annotation_scale(location = "br", width_hint = 0.2) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) 
outmap



#let's check the position of the hits
dftmp <- merge(dftmp, AEZ.lab, by.x="AEZ33", by.y="AEZ")

dapcAGZ <-ggplot(data=dftmp, aes(x=AEZdef, fill=as.factor(DAPC))) + 
  geom_bar(stat="count") +
  scale_x_discrete(n.dod) + 
  scale_fill_d3() +
  labs(x="Agroecological Zone", y="count") +
  guides(fill=guide_legend(title="DAPC \ncluster")) +
  theme_light() +
  theme(legend.position="none")

dapcAGZ

ggsave(outmap, file="output/cowpea.AGZs.clust.jpg", height=6.8, width =5)


save.image(file = "output/metadata.pass.Rdata")


# Prepare admistrative shapes

# get shapes of boardering countries
moz <- getData("GADM", country="MOZ", level=0)
zaf <- getData("GADM", country="ZAF", level=0)
swz  <- getData("GADM", country="SWZ", level=0)
zwe  <- getData("GADM", country="ZWE", level=0)
mwi <- getData("GADM", country="MWI", level=0)
tza <- getData("GADM", country="TZA", level=0)
zmb <- getData("GADM", country="ZMB", level=0)
cod <- getData("GADM", country="COD", level=0)
bwa <-getData("GADM", country="BWA", level=0)
bdi <-getData("GADM", country="BDI", level=0)
rwa <-getData("GADM", country="RWA", level=0)
ken <-getData("GADM", country="KEN", level=0)

#merge all polygons and crop
all <-  rbind (moz, zaf, swz, zwe, mwi, tza, zmb, cod, bwa, bdi, rwa, ken, makeUniqueIDs = TRUE)
all <- crop(all, e)
plot(all)

all_ok <- st_as_sf(all)

save(all, all_ok, file = "output/adm.crop.Rdata")
