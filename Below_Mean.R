##-------------------------------------------
##                  Filtering
##-------------------------------------------
## ----------------------------------------
## option 2: Filtering those below the Mean
## ----------------------------------------
## We are gonna continue with the Option: `Filtering those below the Mean`
## further analysis of the Option 1: `Filtering those below the Mean` is in the script no_filter.R

## Filter those below the mean (0.4891) using 0.4 instead of the complete number
below_MEAN_project_vit_D <- below_MEAN_project_vit_D[, below_MEAN_project_vit_D$assigned_gene_prop > 0.4]
## Save the average expression of a gene in all samples
gene_means_BM <- rowMeans(assay(below_MEAN_project_vit_D, "counts"))
## Delete genes with low mean expression
below_MEAN_project_vit_D <- below_MEAN_project_vit_D[gene_means_BM > 0.1, ]
## Percent of remaining genes:
round(nrow(below_MEAN_project_vit_D) / nrow(unfiltered_vitamin_D) * 100, 2)

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

## ----------------------------------------
##    Differential expression analysis
## ----------------------------------------
## ----------------------------
##            Graph
## ----------------------------
## -------------boxplot---------------
ggplot(as.data.frame(colData(below_MEAN_project_vit_D)),
       aes(x = sra_attribute.treatment,y = assigned_gene_prop, color = sra_attribute.treatment)) +
  geom_boxplot() +
  facet_wrap(~ sra_attribute.cell_type) + # subgraph for each celltype
  theme_bw() +
  labs(title = "BM Distribution of Assigned Gene Prop according to Treatment and Cell Type",
       x = "Treatment", y = "Assigned Gene Prop") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # labels and text

## ----------------------------
## Creating the statistic model
## ----------------------------
## Model as a Function of the variables treatment, cell_type and assigned_gene_prop
models_BM<- model.matrix(~ 0 + sra_attribute.cell_type + sra_attribute.treatment + assigned_gene_prop,
                           data = colData(below_MEAN_project_vit_D)
)
## names of the columns in the matrix
# colnames(models_BM)

## ----------------------------
##            Graph
## ----------------------------
## Dispersion graph of mean-variance trend
vGene_BM <- voom(dge_BM, models_BM, plot = TRUE)
## ----DE-----
## Fitting a linear Model and using an empirical Bayes Method
eb_results_BM <- eBayes(lmFit(vGene_BM))
## Extract a table of the DE genes
de_results_BM <- topTable(
  eb_results_BM,
  coef = 2,
  number = nrow(below_MEAN_project_vit_D),
  sort.by = "none"
)
## Dimensions of the results
## 16
dim(de_results_BM)
##  Count the number of genes with a p-value less than 0.05
table(de_results_BM$adj.P.Val < 0.05)

## ----------------------------------------
##               Results
## ----------------------------------------
## ---------------------------
##            Graph
## ---------------------------
## Relationship between average and the change in gene expression
plotMA(eb_results_BM, coef = 2)
## ----------------------------
##            Graph
## ---------------------------
## ----Volcano Plot----
## Statistical significance and magnitude of the change in expression
## highlight 3 of the most significant genes
volcanoplot(eb_results_BM, coef = 2, highlight = 3, names = de_results_BM$gene_name)

## ----------------------------------------
##          Visualizing DE genes
## ----------------------------------------
## ----------------------------
##            Graph
## ---------------------------
## ----Heatmap----
DE_below_MEan_df <- as.data.frame(colData(
  below_MEAN_project_vit_D)[, c("sra_attribute.cell_type", "sra_attribute.treatment", "assigned_gene_prop")])

colnames(DE_below_MEan_df) <- c("cell_type", "Treatment", "gene_prop")
## get the 50 most significant genes
BM_exprs_heatmap <- vGene_BM$E[rank(de_results_BM$adj.P.Val) <= 50, ]

# identical(rowRanges(below_MEAN_project_vit_D)$gene_id, de_results_BM$gene_id)

## Change IDs to names
gene_names_BM <- rownames(BM_exprs_heatmap)
rownames(BM_exprs_heatmap) <- rowRanges(below_MEAN_project_vit_D)$gene_name[
  match(rownames(BM_exprs_heatmap), rowRanges(below_MEAN_project_vit_D)$gene_id)
]
## Create Heatmap
pheatmap(
  BM_exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = DE_below_MEan_df ## labels
)
