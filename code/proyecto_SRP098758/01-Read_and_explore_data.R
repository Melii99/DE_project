
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


## Usamos el proyecto SRP098758
proj_info <- subset(
  human_projects,
  project == "SRP098758" & project_type == "data_sources"
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

# age - sra_atribute.age
# code - sra_atribute.code
# gender - sra_atribute.gender
# group - sra_attribute.group
# site - sra_attribute.site
# source_name - sra_attribute.source_name
# subjectid - sra_attribute.subjectid
# tissue - sra_attribute.tissue

#################################### iSEE ######################################

library("iSEE")

rowData(rse_gene)
colData(rse_gene)

iSEE::iSEE(rse_gene)

############################### Explorar datos #################################

## Exploramos el atributo de edad
table(rse_gene$sra_attribute.age)

## Exploramos el atributo de code
table(rse_gene$sra_attribute.code)

## Exploramos el atributo de genero
table(rse_gene$sra_attribute.gender)

## Exploramos el atributo de grupo
table(rse_gene$sra_attribute.group)

## Exploramos el atributo de site
table(rse_gene$sra_attribute.site)

## Exploramos el atributo de source_name
table(rse_gene$sra_attribute.source_name)

## Exploramos el atributo de subjectid
table(rse_gene$sra_attribute.subjectid)

## Exploramos el atributo de tejido
table(rse_gene$sra_attribute.tissue)


##### Editamos los datos #####

### No podemos quitar los NAs !!!

### En caso de reemplazar los caracteres por NAs y realizar un droplevels, al intentar
### realizar el an{alisis de expresion diferencial con limma voom() obtenemos un error
### pues las dimensiones de las matrices no coinciden (6 datos no clasif.)

### En caso de retener los NAs para que las dimensiones de las matrices coincidan
### obtenemos un error pues voom() no puede leer los NAs.


rse_gene$sra_attribute.age[rse_gene$sra_attribute.age == "NA"] <- NA
rse_gene$sra_attribute.age <- factor(rse_gene$sra_attribute.age)
rse_gene$sra_attribute.age <- droplevels(rse_gene$sra_attribute.age)

rse_gene$sra_attribute.code[rse_gene$sra_attribute.code == "NA"] <- NA
rse_gene$sra_attribute.code <- factor(rse_gene$sra_attribute.code)
rse_gene$sra_attribute.code <- droplevels(rse_gene$sra_attribute.code)

rse_gene$sra_attribute.gender[rse_gene$sra_attribute.gender == "NA"] <- NA
rse_gene$sra_attribute.gender <- factor(rse_gene$sra_attribute.gender)
rse_gene$sra_attribute.gender <- droplevels(rse_gene$sra_attribute.gender)

rse_gene$sra_attribute.group[rse_gene$sra_attribute.group == "NA"] <- NA
rse_gene$sra_attribute.group <- factor(rse_gene$sra_attribute.group)
rse_gene$sra_attribute.group <- droplevels(rse_gene$sra_attribute.group)

rse_gene$sra_attribute.site[rse_gene$sra_attribute.site == "NA"] <- NA
rse_gene$sra_attribute.site <- factor(rse_gene$sra_attribute.site)
rse_gene$sra_attribute.site <- droplevels(rse_gene$sra_attribute.site)

rse_gene$sra_attribute.source_name[rse_gene$sra_attribute.source_name == "NA"] <- NA
rse_gene$sra_attribute.source_name <- factor(rse_gene$sra_attribute.source_name)
rse_gene$sra_attribute.source_name <- droplevels(rse_gene$sra_attribute.source_name)

rse_gene$sra_attribute.subjectid[rse_gene$sra_attribute.subjectid == "NA"] <- NA
rse_gene$sra_attribute.subjectid <- factor(rse_gene$sra_attribute.subjectid)
rse_gene$sra_attribute.subjectid <- droplevels(rse_gene$sra_attribute.subjectid)

rse_gene$sra_attribute.tissue[rse_gene$sra_attribute.tissue == "NA"] <- NA
rse_gene$sra_attribute.tissue <- factor(rse_gene$sra_attribute.tissue)
rse_gene$sra_attribute.tissue <- droplevels(rse_gene$sra_attribute.tissue)

## Pasar de character a numeric o factor
rse_gene$sra_attribute.age <- as.numeric(rse_gene$sra_attribute.age)
rse_gene$sra_attribute.code <- as.numeric(rse_gene$sra_attribute.code)
rse_gene$sra_attribute.gender <- factor(rse_gene$sra_attribute.gender)
rse_gene$sra_attribute.group <- factor(rse_gene$sra_attribute.group)
rse_gene$sra_attribute.site <- factor(rse_gene$sra_attribute.site)
rse_gene$sra_attribute.source_name <- factor(rse_gene$sra_attribute.source_name)
rse_gene$sra_attribute.subjectid <- factor(rse_gene$sra_attribute.subjectid)
rse_gene$sra_attribute.tissue <- factor(rse_gene$sra_attribute.tissue)

## Resumen de las variables de interés (group, gender, age, site)
summary(as.data.frame(colData(rse_gene)[
  ,
  grepl("^sra_attribute.[group|gender|age|site]", colnames(colData(rse_gene)))
]))

## Factores del atributo grupo: Casos y controles
rse_gene$case <- factor(ifelse(rse_gene$sra_attribute.group == "case (TB)", "case", "control"))
table(rse_gene$case)

## Factores del atributo gender: female y male
rse_gene$sex <- factor(ifelse(rse_gene$sra_attribute.gender == "female", "F", "M"))
table(rse_gene$sex)


## Proporcion de genes para cada grupo (control y caso)
rse_gene$assigned_gene_prop <- rse_gene$recount_qc.gene_fc_count_all.assigned / rse_gene$recount_qc.gene_fc_count_all.total

summary(rse_gene$assigned_gene_prop)

## Diferencias entre los grupos de casos y controles
with(colData(rse_gene), tapply(assigned_gene_prop, case, summary))
