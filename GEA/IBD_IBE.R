library(vegan)
library(geosphere)
setwd("F:/R/IBD/RESULT")
env=read.table("clim.biovalue.txt",head=T)
fst=read.table("fst_matrix.txt",head=T)

### env distance
ENV=as.data.frame(env[,c(4:13)])
ENV<- scale(ENV, center=TRUE, scale=TRUE)
rownames(ENV)=env$POP
env_dist=vegdist(ENV,method="euclidean",binary=FALSE, diag=FALSE, upper=FALSE,na.rm = FALSE)
write.csv(as.matrix(env_dist),file="all_env.csv",quote=F)

#### geographical distance
dist=as.data.frame(env[,c(2:3)])
rownames(dist)=env$POP
muer.dists = distm(dist, fun=distVincentyEllipsoid)
rownames(muer.dists) =env$POP
colnames(muer.dists) =env$POP
write.csv(as.matrix(muer.dists),file="all_dist.csv",quote=F)

#### geo-env dist pearson
env_dist=as.matrix(env_dist)
mantel(muer.dists,env_dist,method="pearson",permutations=999)

####IBD
mantel(fst,muer.dists,method="pearson",permutations=999)

##### IBE
mantel(fst,env_dist,method="pearson",permutations=999)

##### IBE partial
mantel.partial(fst,env_dist,muer.dists,method="pearson",permutations=999)

##### IBD partial
mantel.partial(fst,muer.dists,env_dist,method="pearson",permutations=999)
