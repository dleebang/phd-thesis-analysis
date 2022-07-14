##Evolutionary model for discrete traits

#select disc. variables
data_model_disc <- data %>% 
  select(structure, emission_pattern, note_types)

#z-scale data
data_model_disc <- scale(data_model_disc)

#build a treedata object and extract data with species names for each variable
data_model_disc <- treedata(tree, data_model_disc) %>% .$data
treedata(tree, data_model_disc)$phy
#ER
#ARD
#SYM

#change coding from 0->1 to 1->2 to run fitDiscrete
data_model_disc[,"structure"] <- ifelse(data_model_disc[,"structure"] == 1, 2, 1)
data_model_disc[,"emission_pattern"] <- ifelse(data_model_disc[,"emission_pattern"] == 1, 2, 1)
data_model_disc[,"note_types"] <- case_when(data_model_disc[,"note_types"] == 2 ~ 3,
                                            data_model_disc[,"note_types"] == 1 ~ 2,
                                            data_model_disc[,"note_types"] == 0 ~ 1)


models <- c("ER","ARD", "SYM")

#fitting the models for each variable and each evol. model
#extract aic
evol_model_disc <- apply(data_model_disc, 2, function(x) {
  sapply(models, function(i) {
    fit <- fitDiscrete(tree, dat = x, model = i, transform = "none")
    
    data.frame(AIC = fit$opt$aic, 
         AICc = fit$opt$aicc,
         lnL = fit$opt$lnL)
    
  })
})

evol_model_disc

