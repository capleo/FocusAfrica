############################################################################################
#Genetic DIVERSITY ANALYSIS of Eragrostis teff
############################################################################################

#PRELIMINARY
wd<-"~/focus/cowpea/"
setwd(wd)
rm(list=ls())
options(stringsAsFactors = F)

#load libraries
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(data.table)

#load Identity By Status Mmatrix and its IDs (Plink 1.90)
mibs<-read.table(file = "input/plink.mibs", sep = " ", header = F)
mibs$V345<-NULL

mIDs<-read.table(file = "input/plink.mibs.id", sep = "\t", header = F)
mIDs$V1<- NULL
colnames(mIDs) <- c("IDs")

##insert IDs
mat <- as.matrix(mibs)
rownames(mat)<-mIDs$"IDs"
colnames(mat)<-mIDs$"IDs"

mat[1:10,1:10]

#PLOT a heatmap
pl <- colorRampPalette(c("#ff3419",
                         "#01529f"))

pdf(paste0("output/IBS.pairwise.pruned.pdf"))
heatmap.2(mat, col = pl(50), scale = "none", trace="none", labRow = FALSE, labCol = FALSE)
dev.off()

# # get order of the dendrogram created using "hclust()"
x<- heatmap.2(mat, col = pl(50), scale = "none", trace="none") #class(x)
order <- colnames(x$carpet)
rm(x)

write.table(order, file = "output/order.kinship.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

############################################################################################
#DIVERSITY ANALYSIS PCA & DAPC
############################################################################################

library(adegenet)
library(car)
library(tidyverse)
library(pegas)
library(ape)
library(seqinr)
library(ggplot2)
library(ade4)
library(factoextra)
library(CMplot)
library(ggsci)
library(genetics)
library(EMMREML)
library(compiler)
library(scatterplot3d)
library(hierfstat)
library(poppr)
library(hierfstat)
library(vcfR)
library(data.table)

# load vcf
myvcf <- read.vcfR ("input/snps.maf005.miss09.missind.0.9.heter.site.0.1.heter.ind.0.1.plk.pruned_100_10_05.vcf", verbose = FALSE)

vcfred <- myvcf[i, j, samples = coord$ID, drop =FALSE]


genind <- vcfR2genind(myvcf)

# unique(colnames(myvcf@gt[,1])[duplicated(colnames(myvcf@gt[,1]))]) #check uniqueness

# View(genind)
genind #436 individuals; 2988
summary(genind@loc.n.all) # all loci are biallelic
# genind@ploidy #diploid
# genind_dataframe <- as.data.frame(genind@tab)
# #summary(genind)
# #View(genind_dataframe)
genind_scale <- scaleGen(genind, NA.method="mean") # Impute missing data using mean. This genotype dataset is used for PCA

## 2.2 Calculating basic population genetics statistics
div <- summary(genind) 
class(div$Hexp)
div$Hobs <- as.array(div$Hobs)
div_het <- cbind(div$Hexp, div$Hobs)
div_het <- as.data.frame(div_het) #Create dataframe displaying observed and expected heterozygosity

plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus", cex=0.1, pch=1, col="blue")
#abline(0,1,col="red") #Create expected_het vs. observed_het plot
bartlett.test(list(div$Hexp, div$Hobs)) #assess homogeneity of variance
t.test(div$Hexp, div$Hobs, paired=TRUE, var.equal=FALSE) #assess significance for differences between observed and expected heterozigosity


# 3 PCA
## 3.1 Diversity analysis using PCA
geno.pca<- dudi.pca(genind_scale,cent=TRUE,scale=TRUE,scannf=FALSE, nf=3)
geno.pca_scores <- data.frame(geno.pca$li)# order scores
geno.pca_scores <- data.frame(scale(geno.pca_scores)) #scale scores
fviz_eig(geno.pca) #visualize scree plot
get_eigenvalue(geno.pca) #Visualize loadings statistics (variance explained etc...)
#write.csv(geno.pca_scores, "~/003_SSSUP/teff2/output/pca/geno.pca_scores.csv", row.names=TRUE) #save scores .csv file
#geno.pca_scores <- read.csv("~/003_SSSUP/teff2/output/pca/geno.pca_scores.csv") # load scores .csv file

geno.pca_scores [1:2,1:3]
#geno.pca_scores.df<-data.frame(geno.pca_scores)

#plot 2 PCAs
ggplot(geno.pca_scores, aes(x= Axis1, y = Axis2)) +
  geom_point(size=1) +
  coord_cartesian() +
  theme_minimal()

ggplot(geno.pca_scores, aes(Axis1, Axis3)) +
  geom_point(size=1) + 
  coord_cartesian() +
  theme_minimal()


## DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS
set.seed(9999)

grp <- find.clusters(genind, 
                     #n.pca = 343,
                     #max.n.clust=50,
                     #n.clust = 10
                     ) #kept 400 PCs ##it seems 
dapc1 <- dapc(genind, grp$grp, n.pca = 343, 
              n.da = 5,
              ) #kept 3 discriminant functions
dapc1


y <- as.data.frame(grp$Kstat)
write.table(y, file="output/BIC.vs.nClust.344.PCs.txt", quote = F, sep = "\t")

scatter(dapc1, cell = 0, pch = 18:23, cstar = 2, mstree = FALSE, lwd = 2, lty = 2, posi.da="topright")

##Plot PCA by DAPC assig
assig <- as.data.frame(grp$grp)
geno <- rownames(assig)
rownames(assig) <- NULL
assig <- cbind(geno,assig)
colnames(assig)<-c("geno", "Cluster")

geno1 <- rownames(geno.pca_scores)
rownames(geno.pca_scores) <- NULL
geno.pca_scores<- cbind(geno1,geno.pca_scores)

geno.pca_scores[1:4,1:3]

geno.pca_scores_assig <- merge(x= assig, y= geno.pca_scores, by.x = "geno", by.y= "geno1" )
geno.pca_scores_assig[1:4,1:5]

## save outputs
save.image(file="output/diversity.analysis.step.1.Rdata")
