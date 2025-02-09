##-------------------------------------------
##                  Filtering
##-------------------------------------------
## ----------------------------------------
##        option 1: Don't Filter
## ----------------------------------------
## We are gonna continue with the Option: `Don't Filter`
## further analysis of the Option 2: `Filtering those below the Mean`
## is in the script Below_Mean.R

## Save the average expression of a gene in all samples
gene_means <- rowMeans(assay(project_vit_D, "counts"))
summary(gene_means)

## Delete genes with low mean expression
project_vit_D <- project_vit_D[gene_means > 0.1, ]
## Percent of remaining genes after filtering:
round(nrow(project_vit_D) / nrow(unfiltered_vitamin_D) * 100, 2) ## 57.68%

# Number of samples:
dim(project_vit_D)
## 36834    12

## ----------------------------------------
##          Data Normalization
## ----------------------------------------
## edgeR object creation
dge<-DGEList(
  counts = assay(project_vit_D, "counts"),
  genes = rowData(project_vit_D)
)
## Normalization with the TMM method (default)
dge <- calcNormFactors(dge)
## ----------------------------
##            Graph
## ----------------------------
## -------------boxplot---------------
ggplot(as.data.frame(colData(project_vit_D)),
       aes(x = sra_attribute.treatment,y = assigned_gene_prop, color = sra_attribute.treatment)) +
  geom_boxplot() +
  facet_wrap(~ sra_attribute.cell_type) + # subgraph for each celltype
  theme_bw() +
  labs(title = "Distribution of Assigned Gene Prop according to Treatment and Cell Type",
       x = "Treatment", y = "Assigned Gene Prop") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # labels and text

## ----------------------------
## Creating the statistic model
## ----------------------------
## ----------------------------
##        MODEL MATRIX
## ----------------------------

df <- as.data.frame(colData(project_vit_D)[, c("sra_attribute.cell_type", "sra_attribute.treatment", "assigned_gene_prop")])
colnames(df)<- c("cell_type", "Treatment", "gene_prop")

M_matrix <- ExploreModelMatrix::VisualizeDesign(
  sampleData = df,
  designFormula = ~ 0 + cell_type + Treatment + cell_type:Treatment,
  textSizeFitted = 2
)

cowplot::plot_grid(plotlist = M_matrix$plotlist)

## Model as a Function of the variables treatment, cell_type and assigned_gene_prop
models <- model.matrix(~ 0 + sra_attribute.cell_type + sra_attribute.treatment +
                         sra_attribute.cell_type:sra_attribute.treatment,
                    data = colData(project_vit_D)
                    )
## names of the columns in the matrix
colnames(models)
## ----------------------------------------
##    Differential expression analysis
## ----------------------------------------
## ----------------------------
##            Graph
## ----------------------------
## -----Dispersion graph of mean-variance trend----
vGene <- voom(dge, models, plot = TRUE)
## ----DE-----
## Fitting a linear Model and using an empirical Bayes Method
eb_results <- eBayes(lmFit(vGene))
## Extract a table of the DE genes
de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(project_vit_D),
  sort.by = "none"
)
## Dimensions of the results
## 16
dim(de_results)
##  Count the number of genes with a p-value less than 0.05
table(de_results$adj.P.Val < 0.05)

## ----------------------------------------
##               Results
## ----------------------------------------
## ---------------------------
##            Graph
## ---------------------------
## Relationship between average and the change in gene expression
plotMA(eb_results, coef = 2)
## ----------------------------
##            Graph
## ---------------------------
## ----Volcano Plot----
## Statistical significance and magnitude of the change in expression
## highlight 3 of the most significant genes
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

## ----------------------------------------
##          Visualizing DE genes
## ----------------------------------------
## ----------------------------
##            Graph
## ---------------------------
## ----Heatmap----
DE_no_filter_df <- as.data.frame(colData(project_vit_D)[, c("sra_attribute.cell_type", "sra_attribute.treatment", "assigned_gene_prop")])

colnames(DE_no_filter_df) <- c("cell_type", "Treatment", "gene_prop")
## get the 50 most significant genes
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

identical(rowRanges(project_vit_D)$gene_id, de_results$gene_id)

## Change IDs to names
gene_names <- rownames(exprs_heatmap)
rownames(exprs_heatmap) <- rowRanges(project_vit_D)$gene_name[
  match(rownames(exprs_heatmap), rowRanges(project_vit_D)$gene_id)
]
## Create Heatmap
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = DE_no_filter_df ## labels
)
