## ---------------------------------
##           R_packages
## ---------------------------------
library("recount3")
library(ExploreModelMatrix)
library("edgeR")
## ---------------------------------
##         Data_Download
## ---------------------------------
## Project: SRP124965
human_projects <- available_projects()

## Create a RangedSummarizedExperiment (RSE) object
project_vit_D <- create_rse(subset(human_projects,
  project == "SRP124965" & project_type == "data_sources"
))

## ----------------------------------------
##             Explore_object
## ----------------------------------------

project_vit_D

## class: RangedSummarizedExperiment
## dim: 63856 12
## metadata(8): time_created recount3_version ... annotation recount3_url
##assays(1): raw_counts
## rownames(63856): ENSG00000278704.1 ENSG00000277400.1 ... ENSG00000182484.15_PAR_Y ENSG00000227159.8_PAR_Y
## rowData names(10): source type ... havana_gene tag
## colnames(12): SRR6290091 SRR6290083 ... SRR6290093 SRR6290094
## colData names(175): rail_id external_id ... recount_pred.curated.cell_line BigWigURL

## ------Change_to_read_counts------
assay(project_vit_D, "counts") <- compute_read_counts(project_vit_D)

##-------Facilitate Metadata access------

# project_vit_D$sra.sample_attributes[1:12] # inspect attributes
# expand_sra_attributes(project_vit_D)

project_vit_D <- expand_sra_attributes(project_vit_D)

colData(project_vit_D)[
  ,
  grepl("^sra_attribute", colnames(colData(project_vit_D)))
]

## ----------------------------------------
##                  Model
## ----------------------------------------

## Change character to factor
project_vit_D$sra_attribute.cell_type <- factor(project_vit_D$sra_attribute.cell_type)
project_vit_D$sra_attribute.treatment <- factor(project_vit_D$sra_attribute.treatment)
project_vit_D$sra_attribute.source_name <- factor(project_vit_D$sra_attribute.source_name)

## variables of interest: cell_type, treatment
summary(as.data.frame(colData(project_vit_D)[
  ,
  grepl("^sra_attribute.[cell_type|treatment]", colnames(colData(project_vit_D)))
  ]))

## Proportion of genes
project_vit_D$assigned_gene_prop<-project_vit_D$recount_qc.gene_fc_count_all.assigned / project_vit_D$recount_qc.gene_fc_count_all.total
summary(project_vit_D$assigned_gene_prop)
## Note: The best sample has around 0.6430 of lectures assigned and 75% of the samples have less than 0.5697 of the lectures assigned
## Quality is not the best ...

## ---------
with(colData(project_vit_D), plot(sra_attribute.treatment,assigned_gene_prop))
with(colData(project_vit_D), tapply(assigned_gene_prop, sra_attribute.treatment, summary))
with(colData(project_vit_D), tapply(assigned_gene_prop, sra_attribute.cell_type, summary))
## Note: This...

## save data
unfiltered_vitamin_D<-project_vit_D
# project_vit_D<-unfiltered_vitamin_D # reverse

##-------------------------------------------
##                  Filtering
##-------------------------------------------
hist(project_vit_D$assigned_gene_prop)
## In this step two options were plausible: Don't Filter and Filtering those below the Mean
## further analysis on both can be found at the no_filter.R and Below_Mean.R scripts
