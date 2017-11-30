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


## The three tables in genomics

1. Let's look at some genomic data to see the three tables in action. We might first need to install the necessary packages: 
```r
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase"))

``` 





