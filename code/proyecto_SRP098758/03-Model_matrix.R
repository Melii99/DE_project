######################## Definir modelo estadistico ############################

library("ggplot2")

## Boxplot para comparar por casos y controles
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = case)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Case Group")

## Boxplot para comparar por sexo
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = sex)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Sex")

## Modelo estadistico con model.matrix
mod <- model.matrix(~ sra_attribute.group + sra_attribute.gender + sra_attribute.age + sra_attribute.site +
                      sra_attribute.code + sra_attribute.subjectid, data = colData(rse_gene))

############################## Analisis Exp. Dif. ##############################

## Análisis de expresión diferencial con limma
#library("limma")
#vGene <- voom(dge, mod, plot = TRUE)
