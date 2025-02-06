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
##        option 1: Don't Filter
## ----------------------------------------
## We are gonna continue with the Option: `Don't Filter`
## further analysis of the Option: `Filtering those below the Mean` is in the script Below_Mean_3.R

## option 2: Don't Filter
gene_means <- rowMeans(assay(project_vit_D, "counts"))
summary(gene_means)
project_vit_D <- project_vit_D[gene_means > 0.1, ]
round(nrow(project_vit_D) / nrow(unfiltered_vitamin_D) * 100, 2) ## 57.68

# Number of samples:
(project_vit_D)
## 36834    12
## ----------------------------------------
##          Data Normalization
## ----------------------------------------

dge<-DGEList(
  counts = assay(project_vit_D, "counts"),
  genes = rowData(project_vit_D)
)
dge <- calcNormFactors(dge)

mod <- model.matrix(~ sra_attribute.cell_type + sra_attribute.treatment + assigned_gene_prop,
                    data = colData(project_vit_D)
)
colnames(mod)



vGene <- voom(dge, mod, plot = TRUE)

eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(project_vit_D),
  sort.by = "none"
)

head(de_results)



