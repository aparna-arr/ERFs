library(MASS)
library(caret)

file<-read.delim("lda_data.txt", header=T)
data<-as.data.frame(file)
data.t <- predict(preProcess(data, c("BoxCox", "center", "scale")), data)
r <- lda(formula = sample ~ .,
	data = data.t,
	prior = c(1,1,1)/3)

r
png("lda.png")
plot(r, col=as.integer(data$sample), pch=20, panel = function(x, y, ...) points(x, y, ...)) 

#as.integer(data$sample)

names<-unique(unlist(data$sample))


legend(x="topright", legend=factor(names), pch=20, col=factor(names)) 
 
dev.off()

# plottin the lda:
# https://stat.ethz.ch/pipermail/r-help/2012-May/313988.html
#> png("lda.png")
#> plot(r, abbrev=TRUE, col=as.integer(data$sample))


