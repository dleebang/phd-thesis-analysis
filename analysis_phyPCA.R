  ultra.tree <- ape::read.tree('bok_ultra.tre')
  ggplot2::theme_set(ggplot2::theme_bw())
  
  #select cont. variables
  data_model <- data %>% select(larynx_size:constrictor_post_pcsa) %>% log() %>% scale()

  #build a treedata object and extract data with species names for each variable
  data_model <- treedata(ultra.tree, data_model) %>% .$data


pca <- phyl.pca(ultra.tree, data_model, method = "BM")
plot(pca)
biplot(pca)
plot(pca$Evec[,'PC1'], pca$Evec[,'PC2'])

PC <- as.data.frame(pca$Evec)
ggplot(aes(x = PC[,"PC1"], y = PC[,"PC2"], color = rownames(PC)), data = PC) +
  geom_point()




## EXPLORE PCA
# note: normal pca will be applied to residuals of anatomical traits vs svl
# in this case, phylogeny is already taken into account since we used pgls

# get residuals from anatomic pgls analysis, this time taking muscle vol and fiber length into acc as well
residuals_an <- residuals

# merge residuals with discrete acoustic traits
data_ac <- data[,c('structure', 'emission_pattern', 'note_types')]

data_ac['structure'] <- ifelse(data_ac['structure'] == 0, 'Tonal', 'Pulsado') %>% as.factor()
data_ac['emission_pattern'] <- ifelse(data_ac['emission_pattern'] == 0, 'Continuo', 'Serie') %>% as.factor()
data_ac['note_types'] <- case_when(data_ac['note_types'] == 0 ~ 'Um', 
                                   data_ac['note_types'] == 1 ~ 'Dois_obrigatorio',
                                   data_ac['note_types'] == 2 ~ 'Dois_facultativo') %>% as.factor()

data_pca <- cbind(residuals_an, data_ac)

PC <- prcomp(data_pca[,1:14])
PC <- as.data.frame(PC$x)

ggplot(aes(PC1, PC2, label = rownames(PC), color = data_pca[,'structure']), data = PC) + geom_label(size = 3)

ggplot(aes(PC1, PC2, label =  rownames(PC), color = data_pca[,'emission_pattern']), data = PC) + geom_label(size = 3)

ggplot(aes(PC1, PC2, label =  rownames(PC), color = data_pca[,'note_types']), data = PC) + geom_label(size = 3)
