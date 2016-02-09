.libPaths()

library(MASS)
library(devtools)
library(ggbiplot)
library(caret)
library(cluster)
library(randomForest)

file<-read.delim("lda_data.txt", header=T, row.names=NULL)
data<-as.data.frame(file)


#data.n<-data
#data.n[data.n == 0] <- NA

#log.dat <- log(data.n[,1:length(data.n) - 1])
#data.n <- na.omit(data.n)

# need odd # of ntree for consistent results
rf<-randomForest(data, importance=TRUE, ntree=100001, mtry=5)

importance(rf)

data.n <- data[,length(data)]

dat.t <- predict(preProcess(data, c("BoxCox", "center", "scale")), data)
dat.t<-dat.t[,1:length(dat.t) - 1]
dat.pca<-prcomp(dat.t)

print(dat.pca)


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

table(groups)

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


pdf("pca.pdf")

plot(dat.pca, type="l")

#ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=data.n, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-5,5), xlim=c(-10, 10))

barplot(100*var.explained, las=2, xlab='PCs', names=seq(1:length(dat.t)), ylab='% Variance Explained', col=colors)

# no kmeans
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=data.n, ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-10,10), xlim=c(-5, 15))

plot(pam(dat.t, pamk.best$nc))
# include kmeans
#ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(km$cluster), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))

plot(hc, labels=FALSE)
#plot(pv.1)
#plot(pv.2)
#plot(pv, labels=FALSE)
# hclust with Ward D2 alg
# original biplot (USE)
ggbiplot(dat.pca, obs.scale=1, var.scale=1, groups=factor(groups), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top') + coord_cartesian(ylim=c(-20,15), xlim=c(-10, 25))

	#	ggbiplot(dat.pca, obs.scale=1, var.scale=1, choices=c(i,j), groups=factor(groups), ellipse=TRUE, circle=TRUE) + scale_color_discrete(name='') + theme(legend.direction = 'horizontal', legend.position = 'top')
#temporary
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


