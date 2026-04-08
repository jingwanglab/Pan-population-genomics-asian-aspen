#! /xxx/R/bin/Rscript --no-save --no-restore
.libPaths("/xxx/R/library")
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA
library(data.table)
setwd("/xxx/RDA/snp_indel_sv_RDA")
gen=fread("input/pkdata3.RDA.geno",header=F)
POP=read.table("input/popID.txt",header=F)
LOCI=read.table("input/ID.txt",header=F)
rownames(gen)=as.character(POP$V1)
colnames(gen)=as.character(LOCI$V1)
env <- read.csv("now.csv",head=T)
env$IND <- as.character(env$IND)
str(env) # Look at the structure of the data frame
env$ID<- as.character(env$ID) # Make individual names characters (not factors)
gen.imp =gen
# Confirm that genotypes and environmental data are in the same order
identical(rownames(gen.imp), env[,1]) 
pdf("pairs.panels.env,pdf")
pairs.panels(env[,4:13], scale=T)
dev.off()
pred <- env[,4:13]
pred.pca <- rda(pred, scale=T)
summary(pred.pca)$cont
pdf("pca.env,pdf")
screeplot(pred.pca, main = "Screeplot: Eigenvalues of Wolf Predictor Variables")
## c.rrelations between the PC axis and predictors:
round(scores(pred.pca, choices=1:10, display="species", scaling=0), digits=3)
dev.off()

#run partialrda
wolf.rda <- rda(gen.imp ~ Bio02+Bio03+Bio05+Bio11+Bio13+Bio15+Bio17+SRQ3+WSQ1+CECSOL+ Condition(PC1), data=pred, scale=T)
wolf.rda
write.table(RsquareAdj(wolf.rda),"RsquareAdj", sep="\t", quote=F, row.names=T)
x1=summary(eigenvals(wolf.rda, model = "constrained"))
write.table(x1,"eigenvals", sep="\t", quote=F, row.names=T)
pdf("wolf.rda.sum.pdf")
screeplot(wolf.rda)
dev.off()

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
# if you needed to be very conservative and only identify those loci under very strong selection (i.e., minimize false positive rates), you could increase the number of standard deviations to 3.5.
cand1 <- outliers(load.rda[,1],3.5)
cand2 <- outliers(load.rda[,2],3.5)
cand3 <- outliers(load.rda[,3],3.5)
cand4 <- outliers(load.rda[,4],3.5)
cand5 <- outliers(load.rda[,5],3.5) 
cand6 <- outliers(load.rda[,6],3.5)
cand7 <- outliers(load.rda[,7],3.5)
cand8 <- outliers(load.rda[,8],3.5)
cand9 <- outliers(load.rda[,9],3.5) 
cand10 <- outliers(load.rda[,10],3.5)

ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4) + length(cand5) + length(cand6) + length(cand7) + length(cand8) + length(cand9) + length(cand10)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6,times=length(cand6)), names(cand6), unname(cand6))
cand7 <- cbind.data.frame(rep(7,times=length(cand7)), names(cand7), unname(cand7))
cand8 <- cbind.data.frame(rep(8,times=length(cand8)), names(cand8), unname(cand8))
cand9 <- cbind.data.frame(rep(9,times=length(cand9)), names(cand9), unname(cand9))
cand10 <- cbind.data.frame(rep(10,times=length(cand10)), names(cand10), unname(cand10))
colnames(cand1) <- colnames(cand2) <- colnames(cand3)<- colnames(cand4) <- colnames(cand5) <- colnames(cand6) <- colnames(cand7) <- colnames(cand8) <- colnames(cand9) <- colnames(cand10)<- c("axis","snp","loading")
cand <- rbind(cand1,cand2,cand3,cand4,cand5,cand6,cand7,cand8,cand9,cand10)
cand <- cand[!duplicated(cand)] 
cand$snp <- as.character(cand$snp)
foo <- matrix(nrow=(ncand), ncol=10) 
colnames(foo) <- c("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17","SRQ3","WSQ1","CECSOL")
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo) 
cand <- cand[!duplicated(cand$snp),]

write.table(cand,file="result/cand.txt",quote=F)
#### Next, take the intersection of the results obtained by RDA, lfmm and gemma




