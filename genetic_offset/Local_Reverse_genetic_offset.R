require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
library(sf)

#Read in population locations & climate data
setwd('/xxx/migrant/')
pops <- read.csv("NOW44.csv")
pops =pops [,c(1,2,3)]
colnames(pops)=c('code','lat','long')
#choose predictors
predNames <- c("Bio02","Bio03","Bio05","Bio11","Bio13","Bio15","Bio17")
present=list.files("/xxx/env/present/",full.names = T)
presClim <- stack(present)
presClim <- presClim[[predNames]]
#Creates pred data for gf (cols = population name, long, lat, climate data)
pred <- data.frame(pop=pops$code,long=pops$long, lat=pops$lat, raster::extract(presClim, pops[,c("long","lat")]),stringsAsFactors=FALSE)


#GF model
#read in maf data
snpsAll<- read.csv("11694_frq.txt", stringsAsFactors=FALSE,sep=' ')
#get snp names
snps <- cbind(pops[,c(1,2,3)],snpsAll)
#merge climate data and maf
snps <- merge(pred, snps, by.x="pop", by.y="code", all.x=TRUE)
#run gradient forest model
gfMod <- gradientForest(data=snps, predictor.vars=predNames,response.vars=colnames(snpsAll),
                        ntree=500, 
                        maxLevel=log2(0.368*nrow(snps)/2), trace=T, 
                        corr.threshold=0.70)
#save(gfMod,'gfMod')
#read in shape file
shp <- shapefile("./shp/pd.shp")
#load future climate data by BCC
future_80=dir('./xxx/future/future_BCC_126_2080/',pattern = '.tif',full.names = T)
futClims <- stack(future_80) #stack future climate layers
futClims <- futClims[[predNames]]
futClimDat <- as.data.frame(futClims, xy=TRUE, na.rm=TRUE)
#transform future climate data with gradient forest model
futClimDatGF <- data.frame(futClimDat[,c("x","y")], predict(gfMod,futClimDat[,predNames])) 
#transform current climate data with gradient forest model
pre_clim=mask(presClim, shp)
popDatGF <- as.data.frame(pre_clim, xy=TRUE, na.rm=TRUE)
popDatGF <- data.frame(popDatGF[,c("x","y")], predict(gfMod, popDatGF[,predNames]))
popDatGF <- split(popDatGF, seq(nrow(popDatGF)))

#Forward offset calculation
cl <- makeCluster(40)
registerDoParallel(cl)
forwardOffsetGF <- foreach(i = 1:length(popDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
  #get the focal population
  onePopGF <- popDatGF[[i]]
  #get destination populations and add gf distance
  combinedDatGF <- futClimDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], futClimDatGF[,predNames]))
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  #choose the pixels with the minimum gfOffse
  #############
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  #calculate the distance to the sites with minimum gfOffse, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dist == min(minCoordsGF$dists)),]
  #if multiple sites have the same gfOffset, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  #get local offset
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y ==coordGF$y),"gfOffset"]
  #get the minimum predicted fst - forward offset in this case
  minValGF <- minCoordsGF$gfOffset
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF, forwardOffset=minValGF, predDist=toGoGF, bearing=bearGF,x2=minPtGF[[1]],y2=minPtGF[[2]])
}
stopCluster(cl)
#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#forwardOffset: forward offset
#predDist: distance to site of forward offset
#bearing: bearing to site of forward offset
#x2/y2: coordinate of site of forward offset
forwardOffsetGF <- do.call(rbind, forwardOffsetGF)
write.csv(forwardOffsetGF,paste0("./future_BCC_126_2080_forwardOffsetGF_bio7.csv"), row.names=FALSE)

#Reverse offset calculation
#Getting all coordinates in the range in current climate and transform using gf
popDatGF <- na.omit(as.data.frame(mask(presClim, shp), xy=TRUE))
popDatGF <- data.frame(popDatGF[,c("x","y")], predict(gfMod, popDatGF[,predNames]))

#Gets climate data from the range in future climate and transform using gf
futClimMask <- mask(futClims, mask=shp)
futClimDat <- as.data.frame(futClimMask, xy=TRUE, na.rm=TRUE)
futClimDatGF <- data.frame(futClimDat[,c("x","y")],predict(gfMod,futClimDat[,predNames]))

#Reverse offset calculation
cl <- makeCluster(40)
registerDoParallel(cl)
reverseOffsetGF <- foreach(i = 1:nrow(futClimDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
  #get the focal population in future climate
  onePopGF <- futClimDatGF[i,]
  #make prediction between focal population and current climate
  combinedDatGF <- popDatGF[,c("x","y")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], popDatGF[,predNames]))
  ##Get metrics for the focal population
  #coordinate of focal population
  coordGF <- onePopGF[,c("x","y")]
  #choose the pixels with the minimum offset
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
  #if multiple sites have the same fst, and same distance, one is randomly chosen
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  #get local offset
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y == coordGF$y),"gfOffset"]
  #get the minimum predicted offset - reverse offset in this case
  minValGF <- minCoordsGF$gfOffset
  #get distance and coordinates of site that minimizes fst
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("x","y")]
  #get bearing to the site that minimizes fst
  bearGF <- bearing(coordGF, minPtGF)
  #write out
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]],local=offsetGF, reverseOffset=minValGF, predDist=toGoGF, bearing=bearGF, x2=minPtGF[[1]],y2=minPtGF[[2]])
}
stopCluster(cl)
#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#reverseOffset: reverse offset
#predDist: distance to site of reverse offset
#bearing: bearing to site of reverse offset
#x2/y2: coordinate of site of reverse offset
reverseOffsetGF <- do.call(rbind, reverseOffsetGF)
write.csv(reverseOffsetGF,paste0("./future_BCC_126_2080_reverseOffsetGF_7bio.csv"), row.names=FALSE)

#Average the results of the three modelsand standardize to 0-1  
#Use localOffset as an example 
d1=fread("future_126_2080_allmean.txt")
means <- sapply(d1$localOffset , mean, na.rm=T)
d1$means<-means
d1 <- na.omit(d1)
d1$means <- (d1$means-min(d1$means))/(max(d1$means)-min(d1$means))
#means: Unique contribution of SVs and normalized from 0 to 1
write.csv.csv(d1,,paste0("./future_126_2080_localOffset_mean.csv"), row.names=FALSE)

