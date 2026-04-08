library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(psych)
library(corrplot)
library(ggcorrplot)
setwd("/xxx/corrplot")
raw_data<-read.csv("/xxx/corrplot/now.csv",header=TRUE)
raw_data=raw_data[,5:64]

#raw_data<-raw_data[!duplicated(raw_data,fromLast=TRUE),]
col3 <- colorRampPalette(c("#2082DD", "white", "#FF3F3F")) 
corrmatrix = cor(raw_data, method = "spearman")
rownames(corrmatrix) = as.character(names(raw_data))
colnames(corrmatrix) = as.character(names(raw_data))
res1 <-cor.mtest(corrmatrix, conf.level= .95)
corrplot.mixed(corr=corrmatrix,lower="number",upper="circle",diag="u",upper.col =col3(20),
               lower.col = col3(20),number.cex=0.9,p.mat= res1$p,sig.level= 0.05,
               bg = "white",
               is.corr = TRUE, outline = FALSE, mar = c(0,0,3,0),
               addCoef.col = NULL, addCoefasPercent = FALSE, 
               order = c("original"),
               rect.col = "black", rect.lwd = 1, tl.cex = 1.2,
               tl.col = "black", tl.offset = 0.4, tl.srt = 90,
               cl.cex = 1.1, cl.ratio = 0.2,tl.pos="lt",
               cl.offset = 0.5 )
title(main="Correlation coefficient of 58 environmental variables",cex.main=2.1)