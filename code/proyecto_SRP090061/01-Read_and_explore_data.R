
################################ recount3 ######################################

## Load recount3 R package
library("recount3")

## Vemos el URL actual
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

## Cambiamos el url
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

## Vemos el URL modificado
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)

##################### objeto RSE a partir de un proyecto #######################

## Revisemos todos los proyectos con datos de humano en recount3
human_projects <- available_projects()


## Usamos el proyecto SRP090061
proj_info <- subset(
  human_projects,
  project == "SRP090061" & project_type == "data_sources"
)
## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene <- create_rse(proj_info)
rse_gene

############################ Explorando el objeto ##############################

## Número de genes y muestras
dim(rse_gene)

## IDs de nuestros genes y muestras
head(dimnames(rse_gene))

## Nombres de tablas de cuentas que tenemos (RPKM, CPM, counts, logcounts...)
assayNames(rse_gene)

## Inicio de nuestra tabla de cuentas
head(assay(rse_gene))

## Información de los genes en un objeto de Bioconductor
rowRanges(rse_gene)

## Tabla con información de los genes
rowData(rse_gene)

## Tabla con información de las muestras
colData(rse_gene)

########################## Adaptando la información ############################

## Convertamos las cuentas por nucleotido a cuentas por lectura usando compute_read_counts().
assay(rse_gene, "counts") <- compute_read_counts(rse_gene)

## Facilitamos la lectura de la información del experimento
rse_gene <- expand_sra_attributes(rse_gene)
colData(rse_gene)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene)))
]

################################  SRA-atributes ################################

# control - sra_attribute.control
# cre_line - sra_attribute.cre_line
# days_in_culture - sra_attribute.days_in_culture
# source_name - sra_attribute.source_name
# barcoded - sra_attribute.viral_barcoded

################################ iSEE imágen ###################################

#library("iSEE")

#rowData(rse_gene)
#colData(rse_gene)

#iSEE::iSEE(rse_gene)

############################### Explorar datos #################################

## Exploramos el atributo de control
table(rse_gene$sra_attribute.control)

## Exploramos el atributo de cre_line
table(rse_gene$sra_attribute.cre_line)

## Exploramos el atributo de days_in_culture
table(rse_gene$sra_attribute.days_in_culture)

## Exploramos el atributo de source_name
table(rse_gene$sra_attribute.source_name)

## Exploramos el atributo de viral_barcoded
table(rse_gene$sra_attribute.viral_barcoded)


## Resumen de las variables de interés (control, cre_line, days_in_culture, source_name y viral_barcoded)
summary(as.data.frame(colData(rse_gene)[
  ,
  grepl("^sra_attribute.[control|cre_line|days_in_culture|source_name|viral_barcoded]", colnames(colData(rse_gene)))
]))

## Factores del atributo grupo: Casos y controles
rse_gene$case <- factor(ifelse(rse_gene$sra_attribute.control == "TRUE", "control", "case"))
table(rse_gene$case)

## Factores del atributo source_namer: viral_barcoded y
rse_gene$source <- factor(ifelse(rse_gene$sra_attribute.source_name == "viral_barcoded", "Viral", "CESC"))
table(rse_gene$source)

## Proporcion de genes para cada grupo
rse_gene$assigned_gene_prop <- rse_gene$recount_qc.gene_fc_count_all.assigned / rse_gene$recount_qc.gene_fc_count_all.total

summary(rse_gene$assigned_gene_prop)

## Diferencias entre los grupos
with(colData(rse_gene), tapply(assigned_gene_prop, case, summary))
