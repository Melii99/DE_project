######################## Definir modelo estadistico ############################

library("ggplot2")

## Generamos algunas graficas para ayudarnos a decidir el modelo estadístico

## Boxplot para comparar por casos y controles
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = case)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Case Group")

## Boxplot para comparar por cre_line
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = sra_attribute.cre_line)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Cre cell_line")

## Boxplot para comparar por days_in_culture
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = sra_attribute.days_in_culture)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Days in culture")

## Boxplot para comparar por source
ggplot(as.data.frame(colData(rse_gene)), aes(y = assigned_gene_prop, x = source)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Source")

## Modelo estadistico con model.matrix
mod <- model.matrix(~ case + sra_attribute.cre_line + source + sra_attribute.source_name + sra_attribute.viral_barcoded, data = colData(rse_gene))
