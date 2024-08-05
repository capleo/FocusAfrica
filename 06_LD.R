# LINKAGE DISEQUILIBRIUM
rm(list=ls())

# server dirs
inputdir<-"/home/l.caproni/Documents/cowpea/input/"
LDdir <-"/home/l.caproni/Documents/cowpea/output/LD"


# # local dirs
# inputdir<-"~/focus/cowpea/input/"
# LDdir <-"~/focus/cowpea/output/LD"

setwd(inputdir)
options(stringsAsFactors = F)

## DONE on th server

library(stringr)
library(dplyr)
library(data.table)
library(genetics)
library(LDheatmap)
library(zoo)
library(RColorBrewer)

#A reduced SNP subset excluding SNPs of minor allele frequency <5% is used
#convert VCF t hmp and filter to MAF 5%
#Import hapmap file
data<-read.delim("snps.maf005.miss09.missind.0.9.heter.site.0.1.heter.ind.0.1.hmp", check.names=FALSE)
dim(data)

#trim hapmap
data<-data[,-5:-11]
rownames(data)<-data[,1]
data<-data[,-1:-2]
dim(data)

#revert to mbp
data[,2]<-data[,2]/1000000

#converd to matrix for faster substitutions
dat1<-as.matrix(data)
# dat2 <- dat1[,-c(1:3)]
# dat2[1:10,1:10]
dim(dat1)
dat1[1:10,1:10]

dat.chr.pos <- dat1[,1:2]
dat1 <- dat1[,3:ncol(dat1)]
#substitue characters

dat1<-sub("A", "A/A", dat1)
dat1<-sub("C", "C/C", dat1)
dat1<-sub("G", "G/G", dat1)
dat1<-sub("T", "T/T", dat1)
dat1<-sub("R", "A/G", dat1)
dat1<-sub("Y", "C/T", dat1)
dat1<-sub("S", "G/C", dat1)
dat1<-sub("W", "A/T", dat1)
dat1<-sub("K", "G/T", dat1)
dat1<-sub("M", "A/C", dat1)
dat1<-sub("N", "<NA>", dat1)

# merge chr and pos to snps in the right format
merged <- cbind (dat.chr.pos, dat1)
dim(merged)
merged[1:10, 1:10]

data<-data.frame(merged) #revert to data.frame
#data$chrom <- str_replace(data$chrom, "H", "") #get rid of "H"

data<-data[order(data$chrom, data$pos),] #make sure everything is ordered by chr and pos
dim(data)

#delete all unmapped or contig data
unique(data$chrom)
data <- data[ !(data$chrom %in% c("MYUN",0)), ]
dim(data)

#data[1:10,1:7]#actually check
data<-t(data) #transposing

gen<-makeGenotypes(data) #Make genotypes.
save(gen, file="LDgenotypes.Rdata")

#remove snp position
snpspos<-data.frame(t(data[1:2,]))
#snpspos$chrom <- sub("H","", snpspos$chrom)
snpspos[,1]<-as.numeric(as.character(snpspos[,1]))
snpspos[,2]<-as.numeric(as.character(snpspos[,2]))

save(snpspos, file="snpspos.Rdata")



#create heatmaps of pairwise LD for all SNPs, iterating per chromosome. Long process to run on a personal computer
setwd(LDdir)

snpchr<-split(snpspos, snpspos[,1])

for (c in 1:length(snpchr)){
  print(c)
  tmp<-gen[,colnames(gen) %in% rownames(snpchr[[c]])] #extract from the data frame the SNP data from each chromosome
  #remove all non-genotype columns for now
  zz<-lapply(tmp, class)
  notgen<-grep("character", zz)
  if(length(notgen)>0){
    tmp<-tmp[,-notgen]
  }
  #remove also from snpchr
  snpchr[[c]]<-snpchr[[c]][rownames(snpchr[[c]]) %in% colnames(tmp),]
  print(paste("SNP removed", length(notgen), "out of", ncol(tmp)))
  pdf(paste("CHR", c, "LDheatmap.pdf", sep = "."), width = 9, height = 7)
  tmp<-LDheatmap(tmp, genetic.distances=snpchr[[c]][,2], distance="physical",color=heat.colors(20))
  dev.off()
  save(tmp, file=paste("CHR", c, "LDheatmap.Rdata", sep = "."))
}

######################
#restart from here
# setwd(LDdir)
options(stringsAsFactors = F)

#print out easier to handle files
load("snpspos.Rdata")

for (i in 1:11){
  print(i)
  load(paste("CHR", i, "LDheatmap.Rdata", sep = "."))
  tmpld<-as.matrix(tmp$LDmatrix)
  
  #keep positions and rearrange the dataframe / create a position data frame
  tmpos<-snpspos[rownames(snpspos) %in% colnames(tmpld),]
  tmpos$marker<-rownames(tmpos)
  colnames(tmpos)<-c("chr","position","marker")
  
  #collapse the LD matrix
  m<-cbind(which(!is.na(tmpld),arr.ind = TRUE),na.omit(as.vector(tmpld)))
  colnames(m)<-c("m1","m2","r2")
  tmpos<-tmpos[order(tmpos[,2]),]
  tmpos$idx<-1:nrow(tmpos)
  
  #add position information. For each pairwise comparison this will yield marker names, positions and r2
  newpos<-merge(m, tmpos, by.x="m1", by.y="idx", all.x=T)
  newpos1<-merge(newpos, tmpos, by.x="m2", by.y="idx", all.x=T)
  #drop chromosome data
  newpos1<-subset(newpos1, select=-c(chr.x,chr.y))
  
  
  dist<-abs(newpos1$position.x-newpos1$position.y) #calculate distance as the difference in Mbp between two markers
  newpos1<-cbind(newpos1,dist) #add distance vector to matrix
  #rearrange columns in matrix
  m<-newpos1[,c("m1","position.x","m2","position.y","r2","dist")]
  #order based on dist
  m<-m[order(m[,"dist"]),]
  save(m, file=paste("collapsed_matrix_LD_chr",i,"Rdata", sep="."))
  
}

#set interpolation function based on Hill and Weir (1988) equation based on Marroni (2011) script available at https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/
Er<-function(C_,d){
  length(d)
  res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2))/(n*(2+C_*d)*(11+C_*d)))
  return(res)
}

#set additional parameters
n=336 #n of individuals
exp <- expression(italic(paste(displaystyle(r^2))))

flist<-list()
datashd<-list()

#obtain fpoints by chromosomes
for (j in 1:11){ # n is number of chr, 7 in this case
  print(j)
  load(paste("collapsed_matrix_LD_chr",j,"Rdata", sep="."))
  m1<-m
  m1<-m1[order(m1$dist),]
  #drop variables that will not be used
  m1<-m1[,c("r2","dist")]
  ld<-m1[,"r2"]
  d<-m1[,"dist"]
  nlm<-nls(ld~Er(C_,d[order(d)]),start=c(C_=0.1),control=nls.control(maxiter=100))
  C_<-summary(nlm)$parameters[1]
  fpoints<-Er(C_,d[order(d)])
  tmp<-data.frame(d,fpoints,j)
  flist[[j]]<-data.frame(d,fpoints,j)
  ldd<-0.25
  xpos<-tmp[which((tmp[,2]-ldd)<=0)[1],1]
  datashd[[j]]<-c(max(tmp[,2]), xpos)
}

fpointsall<-do.call(rbind, flist)
# fpointsall[,3]<-as.factor(fpointsall[,3])
datashd<-do.call(rbind, datashd)

colnames(fpointsall)<-c("distance","r2","Chromosome")

fpointsall$distanceKb<-fpointsall$distance*1000

#Export maximum LD and LD decay distance per chromosome to .csv file
colnames(datashd)<-c("max", "Mb")
datashd<- as.data.frame(datashd)
datashd$Kb = datashd$Mb * 1000
write.csv(datashd,"Table.LD.csv")
datashd[11,1:3]
#now upgrade the basic plots using sliding windows
