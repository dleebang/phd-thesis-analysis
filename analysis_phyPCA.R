#select cont. variables
data_model <- data %>% select(call_duration, call_rate, call_freq_pk, call_freq_bw, mean_svl)
#z-scale data
data_model <- scale(data_model)
#build a treedata object and extract data with species names for each variable
data_model <- treedata(tree, data_model) %>% .$data

data_model <- data_model[-28,]
pca <- prcomp(data_model)

plot(pca$x[,"PC1"], pca$x[,"PC2"])

PC <- as.data.frame(pca$x)
ggplot(aes(x = PC[,"PC1"], y = PC[,"PC2"], color = rownames(PC)), data = PC) +
  geom_point()
