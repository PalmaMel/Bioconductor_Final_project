## ---------------------------------
##           R_packages
## ---------------------------------
library("recount3")

## ------------------------------
##         Data_Download
## ------------------------------
## Project: SRP124965
human_projects <- available_projects()

## Create a RangedSummarizedExperiment (RSE) object
project_vit_D <- create_rse(subset(human_projects,
  project == "SRP124965" & project_type == "data_sources"
))

## ------------------------------
##       Explore_object
## ------------------------------

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

## ------------------------
##         Models
## ------------------------

## Change character to factor
project_vit_D$sra_attribute.cell_type <- factor(project_vit_D$sra_attribute.cell_type)
project_vit_D$sra_attribute.treatment <- factor(project_vit_D$sra_attribute.treatment)
project_vit_D$sra_attribute.source_name <- factor(project_vit_D$sra_attribute.source_name)

## variables of interest: cell_type, treatment
summary(as.data.frame(colData(project_vit_D)[
  ,
  grepl("^sra_attribute.[cell_type|treatment]", colnames(colData(project_vit_D)))
  ]))

#-----------------------




