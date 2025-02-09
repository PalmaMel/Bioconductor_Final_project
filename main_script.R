## ---------------------------------
##            Path
## ---------------------------------
setwd("~/Desktop/Bioconductor_Final_project")
## ---------------------------------
##           R_packages
## ---------------------------------
library("recount3")
library("edgeR")
library("limma")
library("ggplot2")
library("pheatmap")
library("ExploreModelMatrix")
## ---------------------------------
##         Data_Download
## ---------------------------------
# Download list of all human experiments
human_projects <- available_projects()

## Create a RangedSummarizedExperiment (RSE) object using the ID
## of the choosen project
project_vit_D <- create_rse(subset(human_projects,
  project == "SRP124965" & project_type == "data_sources"
))

## ----------------------------------------
##             Explore_object
## ----------------------------------------

project_vit_D

## ------Change_to_read_counts------
assay(project_vit_D, "counts") <- compute_read_counts(project_vit_D)

##-------Facilitate Metadata access------
# project_vit_D$sra.sample_attributes[1:12] # inspect attributes, no issues found
expand_sra_attributes(project_vit_D)
## The metadata in the single the SRA (Sequence Read Archive) column is
## separated into several ones in colData
project_vit_D <- expand_sra_attributes(project_vit_D)

## Showcases resultant columns
colData(project_vit_D)[
  ,
  grepl("^sra_attribute", colnames(colData(project_vit_D)))
]

## Change character to factor
project_vit_D$sra_attribute.cell_type <- factor(project_vit_D$sra_attribute.cell_type)
project_vit_D$sra_attribute.treatment <- factor(project_vit_D$sra_attribute.treatment)
project_vit_D$sra_attribute.source_name <- factor(project_vit_D$sra_attribute.source_name)

## Proportion of lectures assigned to genes
project_vit_D$assigned_gene_prop<-project_vit_D$recount_qc.gene_fc_count_all.assigned / project_vit_D$recount_qc.gene_fc_count_all.total
summary(project_vit_D$assigned_gene_prop)
## Note: The best sample has around 0.6430 of lectures assigned and 75% of the samples have less than 0.5697 of the lectures assigned

## save data
unfiltered_vitamin_D<-project_vit_D
## project_vit_D<-unfiltered_vitamin_D # reverse

##-------------------------------------------
##                  Filtering
##-------------------------------------------
hist(project_vit_D$assigned_gene_prop)
## In this step two options were plausible: Don't Filter and Filtering those below the Mean
## further analysis on both can be found at the no_filter.R and Below_Mean.R scripts

## get the mean
summary(unfiltered_vitamin_D$assigned_gene_prop)  ## mean:  0.4891

## Create object for this alternative analysis
below_MEAN_project_vit_D<-project_vit_D


