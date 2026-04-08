library(graph3d)
library(scatterplot3d)
data=read.csv("pd_out.csv",header=TRUE)
color=c('#4dbbd5','#fdbf6f','#4daf4a','#ff7f00','#984ea3','#698ed0')
shape=c(0,5,2,1,3)
data$groups <- as.factor(data$groups)
colors=color[as.numeric(data$groups)]
shapes=shape[as.numeric(data$groups)]
scatterplot3d(data[,2:4],color=colors,pch = shapes,angle=20,main="PCA"))
legend(-2,11, legend=levels(data$groups),col = color,bg = "white",box.lwd = 1,text.font = 0.7, horiz=T,pch=shape)

library(ggplot2)
p1=ggplot(data, aes(x=X43.666, y=X9.727,color=groups,shape=groups))  + geom_point(size=6) +
  scale_color_manual(values=c('#4dbbd5','#fdbf6f','#4daf4a','#698ed0'))+
  scale_shape_manual(values = c(0,5,2,1,3))+
    theme_test()+
    theme(legend.title=element_blank())+
  guides(fill='none')

ggsave(p1,filename='PCA1PCA2.pdf',width = 8,height = 7)
