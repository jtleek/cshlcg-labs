# CSHL Computational Genomics Statistics Lab II

## Get set up

1. Load the data
```r
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()
```
2. Make plots pretty
```r
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)
trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)
```
3. Process the data to make it deal with common problems
```r
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 5,]
```

## Clustering

1. Find the distance between the samples
```r
dist1 = dist(t(edata))
```
2. Look at the distance matrix
```r
install.packages("pheatmap")
library(pheatmap)
pheatmap(as.matrix(dist1),cluster_cols=FALSE,cluster_rows=FALSE)
```
3. Perform a hierarchical clustering
```r
hclust1 = hclust(dist1)
```
4. Plot the clustering ([here are ways to make this pretty](http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning))
```r
plot(hclust1)
```
5. One "prettier" plot - with 4 clusters
```r
install.packages("dendextend")
library(dendextend)
dend = as.dendrogram(hclust1)
dend = color_labels(hclust1,4,col=1:4)
plot(dend)
```
6. Another one with defined groups
```r
labels_colors(dend) = c(rep(1,10),rep(2,9))
plot(dend)
```
7. Do a kmeans clustering
```r
kmeans1 = kmeans(edata,centers=3)
names(kmeans1)
```
8. Plot the cluster means
```r
matplot(t(kmeans1$centers),col=1:3,type="l",lwd=3)
```
9. How many genes in each cluster?
```r
table(kmeans1$cluster)
```
10. Heatmap kmeans clustered
```r
pheatmap(as.matrix(edata)[order(kmeans1$cluster),],
  cluster_cols=F,
  cluster_rows=F)
```
11. This is not deterministic!!!
```r
## You may have to run this a few times
kmeans2 = kmeans(edata,centers=3)
table(kmeans1$cluster,kmeans2$cluster)
```

## Many linear regressions
1. Load some packages we will need
```r
source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
```
2. Load a data set
```r
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]
```
3. Set up the model matrix
```r
mod = model.matrix(~ pdata$strain)
fit = lm.fit(mod,t(edata))
names(fit)
```
4. Look at the fitted coefficients
```r
fit$coefficients[,1]
tidy(lm(as.numeric(edata[1, ]) ~ pdata$strain))
```
5. Plot the fitted coefficients
```r
par(mfrow=c(1,2))
hist(fit$coefficients[1,],breaks=100,col=2,xlab="Intercept")
hist(fit$coefficients[2,],breaks=100,col=2,xlab="Strain")
abline(v=0,lwd=3,col=1)
```
6. Plot the residuals
```r
par(mfrow=c(1,2))
plot(fit$residuals[,1],col=2)
plot(fit$residuals[,2],col=2)
```
7. Create an adjusted model
```r
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_adj = lm.fit(mod_adj,t(edata))
fit_adj$coefficients[,1]
```
8. Fit the model with limma
```r
fit_limma = lmFit(edata,mod_adj)
names(fit_limma)
ebayes_limma = eBayes(fit_limma)
```
9. Compare model fits
```r
plot(ebayes_limma$t[,2],fit_adj$coefficients[2,],col=4,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),col="darkgrey",lwd=3)
```
9. Calculate p-values
```r
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1])$P.Value
```
10. Correct for multiple testing
```r
limma_pvals_adj = topTable(ebayes_limma,number=dim(edata)[1])$adj.P.Val
quantile(limma_pvals)
quantile(limma_pvals_adj)
```
11. Correct for multiple testing with q-value
```r
biocLite("qvalue")
library(qvalue)
qval_limma = qvalue(limma_pvals)
summary(qval_limma)
qval_limma$pi0
```

