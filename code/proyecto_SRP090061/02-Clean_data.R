############################## Limpieza de los datos ###########################

## Guardamos el objeto original
rse_gene_nofilter <- rse_gene

## Creamos un histograma con la proporcion de genes
hist(rse_gene$assigned_gene_prop)

## Eliminamos las muestras malas
table(rse_gene$assigned_gene_prop < 0.3)
rse_gene <- rse_gene[, rse_gene$assigned_gene_prop > 0.3]


## Niveles medios de expresión de los genes en las muestras
gene_means <- rowMeans(assay(rse_gene, "counts"))
summary(gene_means)

## Eliminamos genes
rse_gene <- rse_gene[gene_means > 0.1, ]

## Dimensiones iniciales
dim(rse_gene_nofilter)
## Dimensiones finales
dim(rse_gene)

## Porcentaje de genes retenidos
round(nrow(rse_gene) / nrow(rse_gene_nofilter) * 100, 2)

## Creamos un histograma con la proporcion de genes despues de la limpieza de datos
hist(rse_gene$assigned_gene_prop)

########################## Normalización de datos ##############################

library("edgeR")

dge <- DGEList(
  counts = assay(rse_gene, "counts"),
  genes = rowData(rse_gene)
)

dge <- calcNormFactors(dge)
