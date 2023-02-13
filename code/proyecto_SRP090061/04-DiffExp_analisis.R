############################## Analisis Exp. Dif. ##############################

## Análisis de expresión diferencial con limma
library("limma")

## voom
vGene <- voom(dge, mod, plot = TRUE)

## eBayes
eb_results <- eBayes(lmFit(vGene))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene),
  sort.by = "none"
)

## Dimensiones de los resultados
dim(de_results)

## Observamos los resultados
head(de_results)

## Tabla con la cantidad de genes diferencialmente expresados (por caso-control) con FDR < 5%
table(de_results$adj.P.Val < 0.05)

## Visualizamos los resultados estadísticos
plotMA(eb_results, coef = 2)

## Generamos un volcano plot con los resultados
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)

## Obtenemos m{as informacion de los 3 genes con mas DE
de_results[de_results$gene_name %in% c("DENND2A", "DYNCH1", "RP11-862"), ]

########################### Visualizando genes DE ##############################

## Extraer valores de los genes de interés (50 mayormente DE)
exprs_heatmap <- vGene$E[rank(de_results$adj.P.Val) <= 50, ]

## Creemos una tabla con información de las muestras
df <- as.data.frame(colData(rse_gene)[, c("case", "sra_attribute.cre_line", "sra_attribute.days_in_culture")])
## con nombres de columnas mas amigables
colnames(df) <- c("CaseGroup", "cre_line", "days_in_culture")

## Creamos un heatmap con los 50 genes con mayor DE
library("pheatmap")

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)
