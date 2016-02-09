source("pca.R")
d<-dist(dat.t)
d
head(d)
d<-dist(as.matrix(dat.t))
head(d)
hc<-hclust(d)
pdf("hclst.pdf")
plot(hc)
dev.off()
less(hc)
hc
groups.3 = cutree(hc,2)
groups.3
head(groups.3)
pdf("hclust.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
dev.off()
groups.3 = cutree(hc,3)
pdf("hclust.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
dev.off()
groups.c
groups.3
hc<-hclust(d, method="single")
groups.3 = cutree(hc,3)
pdf("hclust.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
dev.off()
pdf("hclust.pdf")
plot(hc, labels=FALSE)
dev.off()
hc<-hclust(d)
pdf("hclust.pdf")
plot(hc, labels=FALSE)
dev.off()
tables(group.3)
table(group.3)
table(groups.3)
groups.3 = cutree(hc,5)
table(groups.3)
groups.3 = cutree(hc,10)
groups.3 = cutree(hc,5)
groups.3 = cutree(hc,10)
table(groups.3)
hc<-hclust(d, method="ward")
hc<-hclust(d, method="ward.D")
pdf("hclust.pdf")
plot(hc, labels=FALSE)
dev.off()
groups.3 = cutree(hc,2)
table(groups.3)
groups.3 = cutree(hc,3)
table(groups.3)
pdf("clust.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
dev.off()
hc<-hclust(d, method="ward.D2")
groups.3 = cutree(hc,3)
pdf("clust.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
dev.off()
table(groups.3)
pdf("clust.pdf")
plot(hc)
plot(hc, labels=FALSE)
dev.off()
pdf("clust.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
dev.off()
groups.3 = cutree(hc,2)
table(groups.3)
pdf("clust.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
dev.off()
q()
source("pca_randomForest2.R")
rf
?rf
?randomForest
library(igraph)
ig<-graph.adjacency(rf$proximity, mode="undirected", weighted=TRUE)
pdf("graph.pdf")
plot(ig)
dev.off()
fc<-fastgreedy.community(ig)
head(fc)
com<-community.to.membership(g, fc$merges, steps= which.max(fc$modularity)-1)
V(g)$color <- com$membership+1
V(g)$color <- com$membership+1
fg
fc
com<-cutat(ig, steps=7)
com<-cutat(fc, steps=7)
V(g)$color <- com$membership+1
com
V(g)$color <- 6
V(ig)$color <- 6
ig$layout <- layout.fruchterman.reingold
V(g)$color <- com+1
V(ig)$color <- com+1
ig$layout <- layout.fruchterman.reingold
pdf("graph.pdf")
plot(ig, vertex.label=NA)
dev.off()
fc$membership
summmary(fc$membership)
summary(fc$membership)
freq(fc$membership)
frequency(fc$membership)
length(fc$membership)
length(fc$membership[fc$membership == 1])
length(fc$membership[fc$membership == 2])
length(fc$membership[fc$membership == 3])
length(fc$membership[fc$membership == 4])
length(fc$membership[fc$membership == 5])
group5<-cbind(data, fc$membership)
head(group5)
clus
clus<-fc$membership
group5<-cbind(data, clus)
head(group5)
group5$clus<-factor(group5$clus)
data.lda<-predict(preProcess(data, c("BoxCox", "center", "scale")), data)
lda<-lda(formula=group5$clus ~ ., data=data.lda, prior=c(1,1,1,1,1)/5)
lda
pdf("test_lda.pdf")
plot(lda)
dev.off()
plot(lda1, col=as.integer(group1data$groups.3), pch=20, panel=function(x,y, ...) points (x,y, ...))
names<-unique(unlist(group1data$groups.3))
legend(x="topright", legend=factor(names), pch=20, col=factor(names))
pdf("test_lda.pdf")
plot(lda, col=as.integer(group5$clus), pch=20, panel=function(x,y, ...) points (x,y, ...))
names<-unique(unlist(group5$clus))
legend(x="topright", legend=factor(names), pch=20, col=factor(names))
dev.off()
plot(lda, col=as.integer(group5$clus), pch=20, panel=function(x,y, ...) points (x,y, ...))
names<-unique(unlist(group5$clus))
legend(x=0.1, y=0.1, legend=factor(names), pch=20, col=factor(names))
dev.off()
pdf("test_lda.pdf")
plot(lda, col=as.integer(group5$clus), pch=20, panel=function(x,y, ...) points (x,y, ...))
legend(x=0.1, y=0.1, legend=factor(names), pch=20, col=factor(names))
dev.off()
pdf("test_lda.pdf")
plot(lda, col=as.integer(group5$clus), pch=20, panel=function(x,y, ...) points (x,y, ...))
legend(x=2, y=2, legend=factor(names), pch=20, col=factor(names))
dev.off()
lda
group5$clus == 1
length(group5$clus[group5$clus == 1])
length(group5$clus[group5$clus == 2])
length(group5$clus[group5$clus == 3])
length(group5$clus[group5$clus == 4])
length(group5$clus[group5$clus == 5])
group3<-group5[,group5$clus!=1 || group5$clus!=2]
data.lda<-predict(preProcess(group3, c("BoxCox", "center", "scale")), group3)
lda<-lda(formula=group3$clus ~ ., data=data.lda, prior=c(1,1,1)/3)
summary(group3$clus)
group3<-group5[,group5$clus != 1 && group5$clus!=2]
data.lda<-predict(preProcess(group3, c("BoxCox", "center", "scale")), group3)
lda<-lda(formula=group3$clus ~ ., data=data.lda, prior=c(1,1,1)/3)
summary(group3$clus)
group3<-group5[group5$clus != 1 && group5$clus!=2,]
summary(group3$clus)
head(group3)
group3<-group5[,group5$clus == 4 || group5$clus == 3 || group5$clus == 5]
head(group3)
summary(group3$clus)
group3[group3$clus == 1] <- NULL
group3[,group3$clus == 1] <- NULL
group3[group3$clus == 1,] <- NULL
group3[group3$clus == 1] <- NULL
group3[,group3$clus == 1] <- NULL
summary(group3)
group3<-group3[!(group3$clus == "1" || group3$clus == "2")]
summary(group3)
group3<-group3[!(group3$clus == 1 || group3$clus == 2)]
summary(group3)
group3<-group3[!(group3$clus == 1 || group3$clus == 2),]
summary(group3)
group3<-group3[!(group3$clus == "1" || group3$clus == "2"),]
summary(group3)
group3<-group3[!(group3$clus == "1" | group3$clus == "2"),]
summary(group3)
data.lda<-predict(preProcess(group3, c("BoxCox", "center", "scale")), group3)
lda<-lda(formula=group3$clus ~ ., data=data.lda, prior=c(1,1,1)/3)
summary(group3$clus)
group3<-group5[!(group5$clus == "1" | group5$clus == "2"),]
summary(group3$clus)
group3$clus<-factor(group3$clus)
summary(group3$clus)
data.lda<-predict(preProcess(group3, c("BoxCox", "center", "scale")), group3)
lda<-lda(formula=group3$clus ~ ., data=data.lda, prior=c(1,1,1)/3)
lda
pdf("test_lda.pdf")
plot(lda, col=as.integer(group3$clus), pch=20, panel=function(x,y, ...) points (x,y, ...))
names<-unique(unlist(group3$clus))
legend(x=2, y=2, legend=factor(names), pch=20, col=factor(names))
legend(x="topright", legend=factor(names), pch=20, col=factor(names))
dev.off()
pdf("test_lda.pdf")
plot(lda, col=as.integer(group3$clus), pch=20, panel=function(x,y, ...) points (x,y, ...))
legend(x="topright", legend=factor(names), pch=20, col=factor(names))
dev.off()
?randomForest
rf
length(data)
pdf("MDS.pdf")
MDSplot(rf, groups3$clus)
MDSplot(rf, group3$clus)
dev.off()
pdf("MDS.pdf")
MDSplot(rf, group5$clus)
dev.off()
dat.pca
pdf("pca.pdf")
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(group5$clus), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') +
coord_cartesian(ylim=c(-20,20), xlim=c(-20,20))
dev.off()
rhistory("network")
?history
savehistory("network")
