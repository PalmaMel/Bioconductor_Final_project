##-------------------------------------------
##                  Filtering
##-------------------------------------------
## ------------------------------
##           R_packages
## ------------------------------
library("limma")
library("ggplot2")
library("pheatmap")
## ----------------------------------------
## option 2:Filtering those below the Mean
## ----------------------------------------
## We are gonna continue with the Option: `Filtering those below the Mean`
## further analysis of the Option 1: `Filtering those below the Mean` is in the script no_filter.R

## Filter those below the mean (0.4891) using 0.4 instead of the complete number
below_MEAN_project_vit_D <- below_MEAN_project_vit_D[, below_MEAN_project_vit_D$assigned_gene_prop > 0.3]

## Save the average expression of a gene in all samples
gene_means_BM <- rowMeans(assay(below_MEAN_project_vit_D, "counts"))
summary(gene_means_BM)

## Delete genes with low mean expression
below_MEAN_project_vit_D <- below_MEAN_project_vit_D[gene_means_BM > 0.1,]
## Percent of remaining genes:
round(nrow(below_MEAN_project_vit_D) / nrow(unfiltered_vitamin_D) * 100, 2) ## 61.45




# Number of samples:
dim(below_MEAN_project_vit_D)
## 39240     9

## ----------------------------------------
##          Data Normalization
## ----------------------------------------
## edgeR object creation
dge_BM<-DGEList(
  counts = assay(below_MEAN_project_vit_D, "counts"),
  genes = rowData(below_MEAN_project_vit_D)
)

## Normalization with the TMM method (default)
dge_BM <- calcNormFactors(dge_BM)
## ----------------------------
##            Graph
## ----------------------------
## This boxplot ...
ggplot(as.data.frame(colData(below_MEAN_project_vit_D)),
       aes(x = sra_attribute.treatment,y = assigned_gene_prop, color = sra_attribute.treatment)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "black") + # Mean in a cross shape
  facet_wrap(~ sra_attribute.cell_type) + # subgraph for each celltype
  theme_bw() +
  labs(title = "BM Distribution of Assigned Gene Prop according to Treatment and Cell Type",
       x = "Treatment", y = "Assigned Gene Prop") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # labels and text

## ----------------------------
## Creating the statistic model
## ----------------------------
## Model as a Function of the variables treatment, cell_type and assigned_gene_prop
models_BM<- model.matrix(~ sra_attribute.cell_type * sra_attribute.treatment + assigned_gene_prop,
                           data = colData(below_MEAN_project_vit_D)
)
## names of the columns in the matrix
colnames(models_BM)
## ----------------------------------------
##    Differential expression analysis
## ----------------------------------------
## ----------------------------
##            Graph
## ----------------------------
## Dispersion graph of mean-variance trend
vGene_BM <- voom(dge_BM, models_BM, plot = TRUE)



