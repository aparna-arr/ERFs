.libPaths()

library(MASS)
library(devtools)
library(ggbiplot)
library(caret)
library(cluster)
library(randomForest)


set.seed(1)
file<-read.delim("lda_data.txt", header=T, row.names=NULL)
data<-as.data.frame(file)

chromHMMfile<-read.delim("chromHMM_to_ERFs_test.bed", header=T, row.names=NULL)
chrHMM<-as.data.frame(chromHMMfile)
chrHMM<-subset(chrHMM, select=-c(chr, start, end, X))

data<-cbind(data,chrHMM)

data.n <- data$sample


#data.n<-data
#data.n[data.n == 0] <- NA

#log.dat <- log(data.n[,1:length(data.n) - 1])
#data.n <- na.omit(data.n)

# need odd # of ntree for consistent results


dat.t <- predict(preProcess(data, c("BoxCox", "center", "scale")), data)
dat.t<-subset(dat.t, select=-c(sample))

#feature removal cor >= 0.6
#source("recursiveVariablePruning.R")

#dat.t<-featureSelectSeed(dat.t, 0.6,1)

#data<-cbind(data,data.n)



dat.pca<-prcomp(dat.t)

print(dat.pca)

# apply the feature selection results to the rf input dataset
#somehow this stopped working
#data<-data[,!(colnames(data) == setdiff(names(data), names(dat.t)))]
fdata<-cbind(data,data.n)

remove<-setdiff(names(data), names(dat.t))
`%ni%` <- Negate(`%in%`)
data<-subset(data,select=names(data) %ni% remove)


#colnames(data)

set.seed(1)
#rf<-randomForest(data, importance=TRUE, ntree=10001, mtry=5)
rf<-randomForest(data, importance=TRUE, ntree=1000001, mtry=5)

importance(rf)

write.table(rf$proximity, "rf_proximity_matrix.txt", sep="\t", row.names=FALSE)
library(igraph)
ig<-graph.adjacency(rf$proximity, mode="undirected", weighted=TRUE)
fc<-fastgreedy.community(ig, E(ig)$weights)
clusnum<-max(fc$membership)
clusnum
com<-cutat(fc, steps=clusnum)
comClus<-factor(fc$membership)

data.network<-cbind(data, fc$membership)
write.table(data.network, "clusters_networks.txt", sep="\t", row.names=FALSE)

library(fpc)
#pamk.best <- pamk(dat.t)

pamk.best<-pamk(rf$proximity)

cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")

#dat.complete<-cbind(dat.t, data.n)
#km<-kmeans(dat.t, centers=pamk.best$nc)

#head(dat.complete)

library(cluster)
#d<-dist(as.matrix(dat.t))
d<-dist(as.matrix(rf$proximity))

hc<-hclust(d, method="ward.D2")
groups<-cutree(hc, pamk.best$nc)

groups.3<-cutree(hc, 3)

table(groups)
table(groups.3)

group0data<-cbind(data, groups)

group0data$groups<-factor(group0data$groups)

#group0data<-subset(group0data, select=-c(sample))

data0.lda<-predict(preProcess(group0data, c("BoxCox", "center", "scale")), group0data)
lda0<-lda(formula=groups ~ ., data=data0.lda, prior=c(1,1)/2)

lda0

#group1.0<-group0data[group0data$groups == 1,]
#group2.0<-group0data[group0data$groups == 2,]
#
#group1.0<-subset(group1.0, select=-c(groups))
#group2.0<-subset(group2.0, select=-c(groups))
#
#dat.t.1.0 <- predict(preProcess(group1.0, c("BoxCox", "center", "scale")), group1.0)
#dat.pca.1.0<-prcomp(dat.t.1.0)
#
#dat.t.2.0 <- predict(preProcess(group2.0, c("BoxCox", "center", "scale")), group2.0)
##dat.t.2.0<-dat.t.2.0[,1:length(dat.t.2.0) - 2]
#dat.pca.2.0<-prcomp(dat.t.2.0)

group1data <- cbind(data, groups.3)

write.table(group1data, "3group_table_rf_clusters.txt", sep="\t", row.names=FALSE)

group1data$groups.3<-factor(group1data$groups.3)

#group1data<-subset(group1data, select=-c(sample))

data1.lda<-predict(preProcess(group1data, c("BoxCox", "center", "scale")), group1data)
lda1<-lda(formula=groups.3 ~ ., data=data1.lda, prior=c(1,1,1)/3)

lda1

#group1.1<-group1data[group1data$groups.3 == 1,]
#group2.1<-group1data[group1data$groups.3 == 2,]
#group3.1<-group1data[group1data$groups.3 == 3,]
#
#group1.1<-subset(group1.1, select=-c(groups.3))
#group2.1<-subset(group2.1, select=-c(groups.3))
#group3.1<-subset(group3.1, select=-c(groups.3))
#
#dat.t.1.1 <- predict(preProcess(group1.1, c("BoxCox", "center", "scale")), group1.1)
##dat.t.1.1<-dat.t.1.1[,1:length(dat.t.1.1) - 2]
#dat.pca.1.1<-prcomp(dat.t.1.1)
#
#dat.t.2.1 <- predict(preProcess(group2.1, c("BoxCox", "center", "scale")), group2.1)
##dat.t.2.1<-dat.t.2.1[,1:length(dat.t.2.1) - 2]
#dat.pca.2.1<-prcomp(dat.t.2.1)
#
#dat.t.3.1 <- predict(preProcess(group3.1, c("BoxCox", "center", "scale")), group3.1)
##dat.t.3.1<-dat.t.3.1[,1:length(dat.t.3.1) - 2]
#dat.pca.3.1<-prcomp(dat.t.3.1)
#

library(parallel)
library(pvclust)


#pv<-pvclust(t(dat.t), method.hclust="ward.D2",parallel=as.integer(50))
#pv.1<-pvclust(dat.t, parallel=as.integer(25))
#pv.2<-pvclust(dat.t, method.hclust="ward.D2",parallel=as.integer(25))
#groups.pv<-cutree(pv, pamk.best$nc)


var.explained = dat.pca$sdev^2 / sum(dat.pca$sdev^2) 

colors<-character(length(var.explained))
curr_sum <- 0
for (i in 1:length(var.explained)) {
	curr_sum <- curr_sum + var.explained[i]
	if (curr_sum < 0.90) {
		colors[i] <- c("red")
	} else {
		colors[i] <- c("grey")
	}
}


pdf("pca2.pdf")

plot(dat.pca, type="l")

#ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=data.n, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-5,5), xlim=c(-10, 10))

barplot(100*var.explained, las=2, xlab='PCs', names=seq(1:length(dat.t)), ylab='% Variance Explained', col=colors)

# no kmeans
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=data.n, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top')+
coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

plot(pam(dat.t, pamk.best$nc))
# include kmeans
#ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(km$cluster), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))

plot(hc, labels=FALSE)
#plot(pv.1)
#plot(pv.2)
#plot(pv, labels=FALSE)
# hclust with Ward D2 alg
# original biplot (USE)
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top')+ 
coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

#ggbiplot(dat.pca.1.0, obs.scale=1, var.scale=1, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') +
#coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

#ggbiplot(dat.pca.2.0, obs.scale=1, var.scale=1, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') +
#coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

plot(lda0)
#plot(lda0, col=as.integer(group0data$groups), pch=20, panel=function(x,y, ...) points (x,y, ...))
#names<-unique(unlist(group0data$groups))
#legend(x="topright", legend=factor(names), pch=20, col=factor(names))

MDSplot(rf, groups)

ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups.3), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + 
coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

#ggbiplot(dat.pca.1.1, obs.scale=1, var.scale=1, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') +
#coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

#ggbiplot(dat.pca.2.1, obs.scale=1, var.scale=1, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') +
#coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

#ggbiplot(dat.pca.3.1, obs.scale=1, var.scale=1, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') +
#coord_cartesian(ylim=c(-20,20), xlim=c(-20,20)) 

plot(lda1, col=as.integer(group1data$groups.3), pch=20, panel=function(x,y, ...) points (x,y, ...))
names<-unique(unlist(group1data$groups.3))
legend(x="topright", legend=factor(names), pch=20, col=factor(names))

MDSplot(rf, groups.3)

ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(comClus), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') +
coord_cartesian(ylim=c(-20,20), xlim=c(-20,20))


MDSplot(rf, comClus)


	#	ggbiplot(dat.pca, obs.scale=1, var.scale=1, choices=c(i,j), groups=factor(groups), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

#pdf("newpdf.pdf")
#for (i in 1:9) {
#	inc<-i+1
#	for (j in inc:10) {
#		plot<-ggbiplot(dat.pca, obs.scale=1, var.scale=1, choices=c(i,j), groups=factor(groups), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))
#		print(plot)
#	}
#}
#ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=data.n, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') 


