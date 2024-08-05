############################################################################################
# Genotype-environment associations: partial RDA
# Author: Leonardo Caproni
############################################################################################

#PRELIMINARY
wd<-"~/focus/cowpea"
setwd(wd)
rm(list=ls())
options(stringsAsFactors = F)

library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(tidyverse)
library(raster)
library(rgdal)
#library(geodata)
#library(rgeobounadries)
library(ggplot2)
library(data.table)
library(BiodiversityR)
library(dismo)
library(vcfR)
library(robust)
library(adegenet)

#load all data
load("metadata.beforeGF.pruned.set.Rdata")
myvcf <- read.vcfR ("input/snps.maf005.miss09.missind.0.9.heter.site.0.1.heter.ind.0.1.plk.pruned_100_10_05.vcf", verbose = FALSE)
meta.gen <- as.data.frame(myvcf@fix)

# filter out
vcfred <- myvcf[samples = rownames(env.gf.mem), drop =TRUE]
genind <- vcfR2genind(vcfred)

#Genetic Structure as PCA
# impute missing with mean value
genind_scale <- scaleGen(genind, NA.method = "mean") # Impute missing data using mean. This genotype dataset is used for PCA

# Diversity analysis using PCA
geno.pca<- dudi.pca(genind_scale,cent=TRUE,scale=TRUE,scannf=FALSE, nf=3)
geno.pca_scores <- data.frame(geno.pca$li)# order scores
pcagen <- data.frame(scale(geno.pca_scores)) #scale scores
# fviz_eig(geno.pca) #visualize scree plot
# get_eigenvalue(geno.pca) #Visualize loadings statistics (variance explained etc...)

names(pcagen)[1] <- "PC1"
names(pcagen)[2] <- "PC2"
names(pcagen)[3] <- "PC3"

pcagen<-pcagen[order(match(rownames(pcagen), rownames(env.gf.mem))), , drop = FALSE]
write.table(pcagen, file="output/PCA.scores.txt", quote = F, row.names = F, sep = '\t')

#save stuff
save.image(file="output/RDA.metadata.1.Rdata")

###############################################
# Redundancy Analysis 
###############################################
df.rda <- cbind(pcagen, env.gf.mem [,1:7])
# summary(df.rda)

#load genotypes in the correct format for RDA
gen <- read.table("input/snp.cowpea.forGF", header = T, row.names = 1)

# start with geno data, check missingness
perc.miss <- sum(is.na(gen))/(nrow(gen)*ncol(gen)) * 100
perc.miss

# simple, fast imputation (if total missiingness is minor than 5)
gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
gen.imp <- as.data.frame(gen.imp)
#sum(is.na(gen.imp)) # check again, it must be zero

## MATCHING STUFF
# need  geno and env data, same order
ord<-rownames(df.rda)

gen.imp<-subset(gen.imp, rownames(gen.imp) %in% ord)
gen.imp<-gen.imp[order(match(rownames(gen.imp), ord)), , drop = TRUE]

### now subset the other dataframes

## Pure neutral population structure model  
RDA_env <- rda(gen.imp ~ bio15 + bio14 + bio12 + bio2 + bio18 + bio5 + bio3 + Condition(PC1 + PC2 + PC3),  df.rda)
RDA_env

RsquareAdj(RDA_env)
screeplot(RDA_env, main="Eigenvalues of constrained axes")

## load Function rdadapt
#### Function to conduct a RDA based genome scan
rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

# use RDA adapt
rdadapt_env<-rdadapt(RDA_env, 2)

## P-values threshold after Bonferroni correction
thres_env <- 0.05/length(rdadapt_env$p.values)
-log10(thres_env)

outliers <- data.frame(Loci = colnames(gen.imp)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], chr = unlist(lapply(strsplit(colnames(gen.imp)[which(rdadapt_env$p.values<thres_env)], split = "_"), function(x) x[1])))


outliers <- outliers[order(outliers$chr, outliers$p.value),]

## List of outlier names
outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$chr)])


## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = row.names(locus_scores), locus_scores)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "Outliers"
#TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "Outliers"))
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

## Biplot of RDA loci and variables scores
RDA.space <-  ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 2.0) +
  scale_color_manual(values = c("gray90", "#6B4596FF")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 3.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  #facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

# ggsave(RDA.space, file="output/RDA.space.outliers.jpg", height = 5.5, width = 6.5)
# ggsave(RDA.space, file="output/RDA.space.outliers.pdf", height = 5.5, width = 6.5)
# 



Outliers <- rep("Neutral", length(colnames(gen.imp)))
Outliers[colnames(gen.imp)%in%outliers$Loci] <- "Outliers"
#Outliers[colnames(gen.imp)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "Outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(gen.imp)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers)