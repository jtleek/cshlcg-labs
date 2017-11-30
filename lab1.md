# CSHL Computational Genomics Statistics Lab I

## Getting started/projects

1. Open Rstudio
2. Create a new project (slides [46-47 here](http://jtleek.com/advdatasci/slides/03-rmarkdown-and-github.html#46) are potentially helpful)
3. Make sure the package `rmarkdown` is installed. To do this run the command `install.packages("rmarkdown")`. 
3. Create a new R markdown file (slides [4-5 here](http://jtleek.com/advdatasci/slides/03-rmarkdown-and-github.html#4) are potential helpful)
5. Copy and paste [this text](https://raw.githubusercontent.com/SISBID/Module1/gh-pages/labs/rmarkdown-lab.Rmd) into your R markdown file and save it. 
6. Brief R markdown tour (slides [6-13 here are potentially useful](http://jtleek.com/advdatasci/slides/03-rmarkdown-and-github.html#6))
7. Do the tasks described in your R markdown document.
8. Discuss [file naming](http://jtleek.com/advdatasci/slides/02-organizing-version-control-slides.html#25)
9. Discuss project organization(slides [9-14 here are potentially useful](http://jtleek.com/advdatasci/slides/02-organizing-version-control-slides.html#9))
10. Discuss version control (slides [15-23 here](http://jtleek.com/advdatasci/slides/02-organizing-version-control-slides.html#15) are potentially useful.)

## The three tables in genomics

1. Let's look at some genomic data to see the three tables in action. We might first need to install the necessary packages: 
```r
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase"))
``` 
2. Now let's get some old RNA-seq data
```r
con=url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
``` 
3. The object `bodymap.eset` has an expression set you can use. Let's take a look at it
```r
bm = bodymap.eset
bm
``` 
4. The genomics data is here: 
```r
exp_data = exprs(bm)
dim(exp_data)
```
5. Let's just look at a bit of the data
```r
head(exp_data,n=5)
exprs[1:2,1:2]

```
6. Now let's look at the phenotype data
```r 
pheno_data = pData(bm)
dim(pheno_data)
head(pheno_data)
```
7. Finally let's look at the feature data
```r
feature_data = fData(bm)
dim(fData(bm))
fData(bm)[1:10,,1]
``` 

## Some simple summaries

1. Make our life a little easier by breaking out the parts of an ExpressionSet (only recommended for interactive analyses)
```r
pdata = pData(bm)
edata = exprs(bm)
fdata = fData(bm)
```
2. Make some tables to look at variables
```r 
table(pdata$gender)
table(pdata$gender,pdata$race)
```
3. Make a summary of the data
```r
summary(edata)
```
4. Check for missing values using this code
```r
# Use option useNA to include NA's in table
table(pdata$age,useNA="ifany")

# is.na checks for NA values
table(is.na(pdata$age))

# Check for other common missing names
sum(pdata$age==" ")

# Check genomic data for NAs
sum(is.na(edata))

# Make the distribution of NA's by genes
gene_na = rowSums(is.na(edata))
table(gene_na)

# Make the distribution of NA's by samples
sample_na = rowSums(is.na(edata))
table(sample_na)
```
5. Install some packages to make plots pretty :)
```r
devtools::install_github('alyssafrazee/RSkittleBrewer')
library(RSkittleBrewer)
trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)
```

## Data transforms

1. Remember that normal data is nice and symmetric
```r
hist(rnorm(1000),col=2)
```
2. Most genomic data is not
```r
hist(edata[,1],col=2,breaks=100)
```
3. One way to address this skew is to use a transformation. One common transformation is the `log`
```r
hist(log(edata[,1]),col=2,breaks=100)
```
4. But if there are zeros you have some problems
```r
min(log(edata))
```
5. People often remove this problem by adding a small number before taking the log
```r
min(log(edata[,1] + 1))
```
6. Now you can see the zeros in the histogram
```r
hist(log(edata[,1] + 1),breaks=100,col=2)
```
7. Another common choice is log base 2, because then differences between values can be interpreted as "fold changes":
```r
hist(log2(edata[,1] + 1),breaks=100,col=2)
```
8. Another common transform is to remove variables that have little data
```r
hist(rowSums(edata==0),col=2)
```
9. We can find and filter them
```r
low_genes = rowMeans(edata) < 5
table(low_genes)
filt_edata = edata[!low_genes,]
dim(filt_edata)
```
10. Our data is starting to look a bit nicer now
```r
hist(log2(filt_edata[,1] + 1),col=2)
```


## Make some plots


1. Make a boxplot, histogram, or density to look at distributions
```
# boxplot
boxplot(log2(edata+1),col=2,range=0)

# histogram
par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2)

# density 
par(mfrow=c(1,2))
plot(density(log2(edata[,1]+1)),col=2)
lines(density(log2(edata[,2]+1)),col=3)
```
2. Compare distributions with a qq-plot
```r
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3)
```
3. An (often better) way to compare distributions an MA-plot. 
```r
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)
```

## "Tidy" genomic analysis
0. Discuss tidyverse vs. Bioconductor (object oriented)
1. Install the biobroom package
```r
source("https://bioconductor.org/biocLite.R")
biocLite("biobroom")
library(biobroom)
```
2. Let's ["tidy"](https://en.wikipedia.org/wiki/Tidy_data) the expression data
```r
tidybm = tidy(bm,addPheno=TRUE)
tidybm
```
3. Now you can use the ["tidyverse"](https://www.tidyverse.org/)
```r
install.packages("ggplot2")
library(ggplot2)

ggplot(tidybm, aes(x=gender, y=log2(value+1))) +
  geom_boxplot() + ggtitle("Boxplot Showing Effect of Gender on Expression") + theme_bw()
```
4. You can use dplyr to filter (much more on this [here](https://docs.google.com/presentation/d/15meI7W3MeF0afEV5ggdqXfwOIlwy5tcFYTJM-VUGHTs/edit?usp=sharing))
```r
tidybm = tidybm %>% group_by(gene) %>%
  mutate(gene_ave = mean(value)) %>% ungroup()
tidybm %>% filter(gene=="ENSG00000000003") %>%
  select(gene,value,gene_ave)
tidybm = tidybm %>% filter(gene_ave  > 10)
```
5. Make a filtered plot
```r
ggplot(tidybm, aes(x=gender, y=log2(value+1))) +
  geom_boxplot() + ggtitle("Boxplot Showing Effect of Gender on Expression") + theme_bw()
```
6. Make a filtered violin plot
```r
ggplot(tidybm, aes(x=gender, y=log2(value+1),color=gender)) +
  geom_violin() + ggtitle("Boxplot Showing Effect of Gender on Expression") + theme_bw()
```