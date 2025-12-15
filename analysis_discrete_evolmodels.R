##Evolutionary model for discrete traits


#build ordered and directional models of evolution for two state traits
ordered.model.2states<-matrix(c(
  0,1,2,0),2,2,byrow=TRUE,
  dimnames=list(0:1,0:1))


directional.model.2states.01<-matrix(c(
  0,1,0,0),2,2,byrow=TRUE,
  dimnames=list(0:1,0:1))


directional.model1.2states.10 <- matrix(c(
  0,0,1,0),2,2,byrow=TRUE,
  dimnames=list(0:1,0:1))



#testing ordered and directional models of evolution for three state traits
ordered.model.3states<-matrix(c(
  0,1,0,
  2,0,3,
  0,4,0),3,3,byrow=TRUE,
  dimnames=list(0:2,0:2))


directional.model.3states.012<-matrix(c(
  0,1,0,
  0,0,2,
  0,0,0),3,3,byrow=TRUE,
  dimnames=list(0:2,0:2))


directional.model.3states.210<-matrix(c(
  0,0,0,
  1,0,0,
  0,2,0),3,3,byrow=TRUE,
  dimnames=list(0:2,0:2))





#### ANATOMIC

write.tree(ultra.tree, 'bok_ultra.tre')
ultra.tree <- ape::read.tree('bok_ultra.tre')

#select disc. variables
data_model_disc <- data %>% 
  select(base_da_aritenoide, apice_da_aritenoide) %>% 
  mutate(base_da_aritenoide  = case_when(data$base_da_aritenoide == "espesso" ~ 1,
                                         data$base_da_aritenoide == "espesso medio" ~ 2,
                                         data$base_da_aritenoide == "espesso com proeminencia" ~ 3)) %>% 
  mutate(apice_da_aritenoide = case_when(data$apice_da_aritenoide == "estreito" ~ 1,
                                         data$apice_da_aritenoide == "medio" ~ 2,
                                         data$apice_da_aritenoide == "ampla" ~ 3))


#build a treedata object and extract data with species names for each variable
data_model_disc <- treedata(ultra.tree, data_model_disc) %>% .$data



#fit models for base_da_aritenoide
fit_ER_ba <- fitDiscrete(ultra.tree, dat = data_model_disc[,'base_da_aritenoide'], model = 'ER')
fit_ARD_ba <- fitDiscrete(ultra.tree, dat = data_model_disc[,'base_da_aritenoide'], model = 'ARD')
fit_SYM_ba <- fitDiscrete(ultra.tree, dat = data_model_disc[,'base_da_aritenoide'], model = 'SYM')
fit_ordered_ba <- fitDiscrete(ultra.tree, dat = data_model_disc[,'base_da_aritenoide'], model = ordered.model.3states, surpressWarnings=TRUE)
fit_directional_ba <- fitDiscrete(ultra.tree, dat = data_model_disc[,'base_da_aritenoide'], model = directional.model.3states.012, surpressWarnings=TRUE)
fit_directional_ba1 <- fitDiscrete(ultra.tree, dat = data_model_disc[,'base_da_aritenoide'], model = directional.model.3states.210, surpressWarnings=TRUE)

round(
  aic.w(c(AIC(fit_ER_ba), AIC(fit_ARD_ba), AIC(fit_SYM_ba), AIC(fit_ordered_ba), AIC(fit_directional_ba), AIC(fit_directional_ba1))),
  3)

library(lmtest)
lrtest(fit_ER_ba, fit_ordered_ba)


#fit models for apice_da_aritenoide
fit_ER_aa <- fitDiscrete(ultra.tree, dat = data_model_disc[,'apice_da_aritenoide'], model = 'ER')
fit_ARD_aa <- fitDiscrete(ultra.tree, dat = data_model_disc[,'apice_da_aritenoide'], model = 'ARD')
fit_SYM_aa <- fitDiscrete(ultra.tree, dat = data_model_disc[,'apice_da_aritenoide'], model = 'SYM')
fit_ordered_aa <- fitDiscrete(ultra.tree, dat = data_model_disc[,'apice_da_aritenoide'], model = ordered.model.3states, surpressWarnings=TRUE)
fit_directional_aa <- fitDiscrete(ultra.tree, dat = data_model_disc[,'apice_da_aritenoide'], model = directional.model.3states.012, surpressWarnings=TRUE)
fit_directional_aa1 <- fitDiscrete(ultra.tree, dat = data_model_disc[,'apice_da_aritenoide'], model = directional.model.3states.210, surpressWarnings=TRUE)

round(
  aic.w(c(AIC(fit_ER_aa), AIC(fit_ARD_aa), AIC(fit_SYM_aa), AIC(fit_ordered_aa), AIC(fit_directional_aa), AIC(fit_directional_aa1))),
  3)















#### ACOUSTIC



#select disc. variables
data_model_disc <- data %>% 
  select(structure, emission_pattern, note_types)

#z-scale data
#data_model_disc <- scale(data_model_disc)

#build a treedata object and extract data with species names for each variable
data_model_disc <- treedata(ultra.tree, data_model_disc) %>% .$data



#change coding from 0->1 to 1->2 to run fitDiscrete
data_model_disc[,"structure"] <- ifelse(data_model_disc[,"structure"] == 1, 2, 1)
data_model_disc[,"emission_pattern"] <- ifelse(data_model_disc[,"emission_pattern"] == 1, 2, 1)
data_model_disc[,"note_types"] <- case_when(data_model_disc[,"note_types"] == 2 ~ 3,
                                            data_model_disc[,"note_types"] == 1 ~ 2,
                                            data_model_disc[,"note_types"] == 0 ~ 1)






#fit models for structure
fit_ER_stru <- fitDiscrete(ultra.tree, dat = data_model_disc[,'structure'], model = 'ER')
fit_ARD_str <- fitDiscrete(ultra.tree, dat = data_model_disc[,'structure'], model = 'ARD')
fit_ordered_stru <- fitDiscrete(ultra.tree, dat = data_model_disc[,'structure'], model = ordered.model, surpressWarnings=TRUE)
fit_directional_stru <- fitDiscrete(ultra.tree, dat = data_model_disc[,'structure'], model = directional.model, surpressWarnings=TRUE)
fit_directional_stru1 <- fitDiscrete(ultra.tree, dat = data_model_disc[,'structure'], model = directional.model1, surpressWarnings=TRUE)

round(
aic.w(c(AIC(fit_ER_stru), AIC(fit_ARD_str), AIC(fit_directional_stru), AIC(fit_directional_stru1))),
3)


#fit models for emission pattern
fit_ER_emis <- fitDiscrete(ultra.tree, dat = data_model_disc[,'emission_pattern'], model = 'ER')
fit_ARD_emis <- fitDiscrete(ultra.tree, dat = data_model_disc[,'emission_pattern'], model = 'ARD')
fit_ordered_emis <- fitDiscrete(ultra.tree, dat = data_model_disc[,'emission_pattern'], model = ordered.model, surpressWarnings=TRUE)
fit_directional_emis <- fitDiscrete(ultra.tree, dat = data_model_disc[,'emission_pattern'], model = directional.model, surpressWarnings=TRUE)
fit_directional_emis1 <- fitDiscrete(ultra.tree, dat = data_model_disc[,'emission_pattern'], model = directional.model1, surpressWarnings=TRUE)


round(
aic.w(c(AIC(fit_ER_emis), AIC(fit_ARD_emis), AIC(fit_directional_emis), AIC(fit_directional_emis1))),
3)




#fit models for note types
fit_ER_note <- fitDiscrete(ultra.tree, dat = data_model_disc[,'note_types'], model = 'ER')
fit_ARD_note <- fitDiscrete(ultra.tree, dat = data_model_disc[,'note_types'], model = 'ARD')
fit_SYM_note <- fitDiscrete(ultra.tree, dat = data_model_disc[,'note_types'], model = 'SYM')
fit_ordered_note <- fitDiscrete(ultra.tree, dat = data_model_disc[,'note_types'], model = ordered.model1, surpressWarnings=TRUE)
fit_directional_note <- fitDiscrete(ultra.tree, dat = data_model_disc[,'note_types'], model = directional.model1, surpressWarnings=TRUE)
fit_directional_note1 <- fitDiscrete(ultra.tree, dat = data_model_disc[,'note_types'], model = directional.model2, surpressWarnings=TRUE)

round(
aic.w(c(AIC(fit_ER_note), AIC(fit_ARD_note), AIC(fit_SYM_note), AIC(fit_ordered_note), AIC(fit_directional_note), AIC(fit_directional_note1))),
3)

library(lmtest)
lrtest(fit_ER_note, fit_ARD_note)



















###### FORMER APPROACH ######
models <- c("ER","ARD", "SYM")

#fitting the models for each variable and each evol. model
#extract aic
evol_model_disc <- apply(data_model_disc, 2, function(x) {
  sapply(models, function(i) {
    fit <- fitDiscrete(ultra.tree, dat = x, model = i, transform = "none")
    
    data.frame(AIC = fit$opt$aic, 
               AICc = fit$opt$aicc,
               lnL = fit$opt$lnL)
    
  })
})


#get aic weights
evol_model_disc <- lapply(evol_model_disc, function(x) {
  x <- rbind(x, round(aic.w(unlist(x[1,])), 3))
})

evol_model_disc
