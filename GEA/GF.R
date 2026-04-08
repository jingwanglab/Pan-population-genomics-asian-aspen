#! /xxx/envs/R/bin/Rscript --no-save --no-restore
.libPaths("/xxx/envs/R/lib/R/library")
library(gradientForest)
library(data.table)
setwd("/xxx/GF")
env=read.table("/xxx/environ_gf.txt", header=T)
envGF=env[, c("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17","SRQ3","WSQ1","CECSOL")]
pops=read.table("/xxx/GF_lon_lat.txt", header=T)
sites=pops[, c("lon", "lat")]
Grid=fread("/xxx/present.txt",header=T)
grid=Grid[, c("lon", "lat")]
greengrid=Grid[,c("lon", "lat","Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17","SRQ3","WSQ1","CECSOL")]
all_SNPs=read.table("MAF.txt",head=T)
preds <- colnames(envGF)
specs <- colnames(all_SNPs)
nSites <- dim(envGF)[1]
nSpecs <- dim(all_SNPs)[2]
maxLevel <- log2(0.368*nrow(envGF)/2)
all_gfmod <- gradientForest(cbind(envGF, all_SNPs), predictor.vars=colnames(envGF),response.vars=colnames(all_SNPs), ntree=500, compact=T, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)
pdf(file="picture/all_predictoroverallimportance.pdf")
plot.gradientForest(all_gfmod,plot.type="O")
dev.off()

write.table(all_gfmod$Y, file="result/all_gf_Y.txt")
write.table(all_gfmod$X, file="result/all_gf_X.txt")
write.table(all_gfmod$imp.rsq, file="result/all_gf_impRsq.txt")
write.table(all_gfmod$result, file="result/all_gf_result.txt")
write.table(all_gfmod$res.u, file="result/all_gf_res_u.txt")
write.table(all_gfmod$res, file="result/all_gf_res.txt")

#R2
pdf(file="picture/all_R2.pdf")
plot(all_gfmod, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)
dev.off()
