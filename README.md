# GeneCoexp
Applies dimension reduction and clustering techniques and creates interactive clustering visualization plots.

# Install

```r
# installing from bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GeneCoexp")

# installing from github
library(devtools)
install_github('BioinfoUSAL/GeneCoexp')
```

# Usage

```r
library(GeneCoexp)

# single
res <- nrcor(iris[,5],t(iris[,1:4]))

# multi
obj <- multinrcor(t(iris[,1:4]))
plot(obj)
```

