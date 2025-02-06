##-------------------------------------------
##                  Filtering
##-------------------------------------------
## ---------------------------------
##           R_packages
## ---------------------------------
library("recount3")
library(ExploreModelMatrix)
library("edgeR")
library("limma")
## ----------------------------------------
## option 2: Filtering those below the Mean
## ----------------------------------------
below_MEAN_project_vit_D<-project_vit_D
below_MEAN_project_vit_D <- below_MEAN_project_vit_D[, below_MEAN_project_vit_D$assigned_gene_prop > 0.4]
gene_means <- rowMeans(assay(below_MEAN_project_vit_D, "counts"))
summary(gene_means)
below_MEAN_project_vit_D <- below_MEAN_project_vit_D[gene_means > 0.1,]
round(nrow(below_MEAN_project_vit_D) / nrow(unfiltered_vitamin_D) * 100, 2) ## 61.45

# Number of samples:
dim(below_MEAN_project_vit_D)
## 39240     9

