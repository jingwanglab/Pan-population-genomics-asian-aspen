#! /data/apps/R/4.0.3/bin/Rscript --no-save --no-restore
#.libPaths("/xxx/R/x86_64-pc-linux-gnu-library/4.0")
library(raster)
library(geosphere)
library(gdm)
library(foreach)
library(parallel)
library(doParallel)
library(gradientForest)
library(fields)
library(sf)
library(data.table)
setwd("/xxx/GF")
offset_outdir="/xxx/offset"
all_SNPs=read.table("/xxx/11694_freq.txt",head=T)
all_SNPs=all_SNPs[,-1]
env=read.table("/usr_storage2/syp/LZQ/env_data/now/PL20_NOW.csv",head=T)
future_dir <- commandArgs(trailingOnly = TRUE)
envGF=env[, c(("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17"))]
predNames=c(("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17"))
sites=env[, c("lon", "lat")]
Grid=fread("/xxx/NOW.csv",header=T)
greengrid=Grid[,c("lon","lat","Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17")]
grid=Grid[, c("lon", "lat")]

pops=env[,1:3]

preds <- colnames(envGF)
specs <- colnames(all_SNPs)
nSites <- dim(envGF)[1]
nSpecs <- dim(all_SNPs)[2]
maxLevel <- log2(0.368*nrow(envGF)/2)
all_gfmod <- gradientForest(cbind(envGF, all_SNPs), predictor.vars=colnames(envGF),response.vars=colnames(all_SNPs), ntree=500, compact=T, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)
shp <- shapefile("./shp/pd.shp")
future_ouya=dir(future_dir,pattern = '.tif',full.names = T)
all_tgrid=cbind(greengrid[,c("lon","lat")], predict(all_gfmod,greengrid[,c("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17")]))
futClims <- stack(future_ouya) 
futClims <- futClims[[predNames]]
a=unlist(strsplit(future_dir,split="/"))
model=a[7]
SSP=a[8]
out_index=paste(offset_outdir,"/",model,"/",model,"_",SSP,"_pd_forward_offset",sep="")

#transform future climate data with gradient forest model
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)
futClimDatGF <- data.frame(futClimDat[,c("x","y")], predict(all_gfmod,futClimDat[,predNames])) 

#transform current climate data with gradient forest model
forwardOffsetGF=vector(mode='list',length=5)
forwardOffsetGF[[1]]=list()
forwardOffsetGF[[2]]=list()
forwardOffsetGF[[3]]=list()
forwardOffsetGF[[4]]=list()
forwardOffsetGF[[5]]=list()

popDatGF=data.frame(all_tgrid)
colnames(popDatGF)=c("x","y","Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17")
for (i in 1:nrow(popDatGF)){
  onePopGF <- popDatGF[i,]
  print(i)
  #get destination populations and add gf distance
  combinedDatGF <- futClimDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,c("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17")], futClimDatGF[,c("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17")]))
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  #choose the pixels with the minimum gfOffse
  combinedDatGF['dists']=distGeo(p1=coordGF, p2=combinedDatGF[,1:2])
  combinedDatGF_100=combinedDatGF[combinedDatGF['dists']<100000,]
  combinedDatGF_250=combinedDatGF[combinedDatGF['dists']<250000,]
  combinedDatGF_500=combinedDatGF[combinedDatGF['dists']<500000,]
  combinedDatGF_1000=combinedDatGF[combinedDatGF['dists']<1000000,]
  combinedDatGF_unlimited=combinedDatGF
  
  minCoordsGF_100 <- combinedDatGF_100[which(combinedDatGF_100$gfOffset == min(combinedDatGF_100$gfOffset)),]
  minCoordsGF_250 <- combinedDatGF_250[which(combinedDatGF_250$gfOffset == min(combinedDatGF_250$gfOffset)),]
  minCoordsGF_500 <- combinedDatGF_500[which(combinedDatGF_500$gfOffset == min(combinedDatGF_500$gfOffset)),]
  minCoordsGF_1000 <- combinedDatGF_1000[which(combinedDatGF_1000$gfOffset == min(combinedDatGF_1000$gfOffset)),]
  minCoordsGF_unlimited <- combinedDatGF_unlimited[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  #calculate the distance to the sites with minimum gfOffse, and selct the one with the shortest distance
  #minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF_100 <- minCoordsGF_100[which(minCoordsGF_100$dist == min(minCoordsGF_100$dists)),]
  minCoordsGF_250 <- minCoordsGF_250[which(minCoordsGF_250$dist == min(minCoordsGF_250$dists)),]
  minCoordsGF_500 <- minCoordsGF_500[which(minCoordsGF_500$dist == min(minCoordsGF_500$dists)),]
  minCoordsGF_1000 <- minCoordsGF_1000[which(minCoordsGF_1000$dist == min(minCoordsGF_1000$dists)),]
  minCoordsGF_unlimited <- minCoordsGF_unlimited[which(minCoordsGF_unlimited$dist == min(minCoordsGF_unlimited$dists)),]
  #if multiple sites have the same gfOffset, and same distance, one is randomly chosen
  minCoordsGF_100 <- minCoordsGF_100[sample(1:nrow(minCoordsGF_100),1),]
  minCoordsGF_250 <- minCoordsGF_250[sample(1:nrow(minCoordsGF_250),1),]
  minCoordsGF_500 <- minCoordsGF_500[sample(1:nrow(minCoordsGF_500),1),]
  minCoordsGF_1000 <- minCoordsGF_1000[sample(1:nrow(minCoordsGF_1000),1),]
  minCoordsGF_unlimited <- minCoordsGF_unlimited[sample(1:nrow(minCoordsGF_unlimited),1),]
  #get local offset
  offsetGF_100 <- combinedDatGF_100[which(combinedDatGF_100$x == coordGF$x & combinedDatGF_100$y ==coordGF$y),"gfOffset"]
  offsetGF_250 <- combinedDatGF_250[which(combinedDatGF_250$x == coordGF$x & combinedDatGF_250$y ==coordGF$y),"gfOffset"]
  offsetGF_500 <- combinedDatGF_500[which(combinedDatGF_500$x == coordGF$x & combinedDatGF_500$y ==coordGF$y),"gfOffset"]
  offsetGF_1000 <- combinedDatGF_1000[which(combinedDatGF_1000$x == coordGF$x & combinedDatGF_1000$y ==coordGF$y),"gfOffset"]
  offsetGF_unlimited <- combinedDatGF_unlimited[which(combinedDatGF_unlimited$x == coordGF$x & combinedDatGF_unlimited$y ==coordGF$y),"gfOffset"]
  #get the minimum predicted fst - forward offset in this case
  minValGF_100 <- minCoordsGF_100$gfOffset
  minValGF_250 <- minCoordsGF_250$gfOffset
  minValGF_500 <- minCoordsGF_500$gfOffset
  minValGF_1000 <- minCoordsGF_1000$gfOffset
  minValGF_unlimited <- minCoordsGF_unlimited$gfOffset
  #get distance and coordinates of site that minimizes fst
  toGoGF_100 <- minCoordsGF_100$dists
  toGoGF_250 <- minCoordsGF_250$dists
  toGoGF_500 <- minCoordsGF_500$dists
  toGoGF_1000 <- minCoordsGF_1000$dists
  toGoGF_unlimited <- minCoordsGF_unlimited$dists
  minPtGF_100 <- minCoordsGF_100[,c("x","y")]
  minPtGF_250 <- minCoordsGF_250[,c("x","y")]
  minPtGF_500 <- minCoordsGF_500[,c("x","y")]
  minPtGF_1000 <- minCoordsGF_1000[,c("x","y")]
  minPtGF_unlimited <- minCoordsGF_unlimited[,c("x","y")]
    #get bearing to the site that minimizes fst
  bearGF_100 <- bearing(coordGF, minPtGF_100)
  bearGF_250 <- bearing(coordGF, minPtGF_250)
  bearGF_500 <- bearing(coordGF, minPtGF_500)
  bearGF_1000 <- bearing(coordGF, minPtGF_1000)
  bearGF_unlimited <- bearing(coordGF, minPtGF_unlimited)
  #write out
  outGF_100 <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF_100, forwardOffset=minValGF_100, predDist=toGoGF_100, bearing=bearGF_100,x2=minPtGF_100[[1]],y2=minPtGF_100[[2]])
  outGF_250 <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF_250, forwardOffset=minValGF_250, predDist=toGoGF_250, bearing=bearGF_250,x2=minPtGF_250[[1]],y2=minPtGF_250[[2]])
  outGF_500 <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF_500, forwardOffset=minValGF_500, predDist=toGoGF_500, bearing=bearGF_500,x2=minPtGF_500[[1]],y2=minPtGF_500[[2]])
  outGF_1000 <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF_1000, forwardOffset=minValGF_1000, predDist=toGoGF_1000, bearing=bearGF_1000,x2=minPtGF_1000[[1]],y2=minPtGF_1000[[2]])
  outGF_unlimited <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF_unlimited, forwardOffset=minValGF_unlimited, predDist=toGoGF_unlimited, bearing=bearGF_unlimited,x2=minPtGF_unlimited[[1]],y2=minPtGF_unlimited[[2]])
  forwardOffsetGF[[1]][[i]]=outGF_100
  forwardOffsetGF[[2]][[i]]=outGF_250
  forwardOffsetGF[[3]][[i]]=outGF_500
  forwardOffsetGF[[4]][[i]]=outGF_1000
  forwardOffsetGF[[5]][[i]]=outGF_unlimited
 }

#stopCluster(cl)
write.csv( do.call(rbind, forwardOffsetGF[[1]]),paste(out_index,"100km.csv",sep="_"), row.names=FALSE)
write.csv( do.call(rbind, forwardOffsetGF[[2]]),paste(out_index,"250km.csv",sep="_"), row.names=FALSE)
write.csv( do.call(rbind, forwardOffsetGF[[3]]),paste(out_index,"500km.csv",sep="_"), row.names=FALSE)
write.csv( do.call(rbind, forwardOffsetGF[[4]]),paste(out_index,"1000km.csv",sep="_"), row.names=FALSE)
write.csv( do.call(rbind, forwardOffsetGF[[5]]),paste(out_index,"unlimited.csv",sep="_"), row.names=FALSE)