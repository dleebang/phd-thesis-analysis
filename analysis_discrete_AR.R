#compute ancestral state reconstructions based on best fit model (ARD, SYM, ER)
#compute ancestral state reconstructions based on all discrete evol models and compare likelihoods?
#https://www.r-phylo.org/wiki/HowTo/InferringModelsForDiscreteData
#
#
#make.simmap (phytools) - stochastic character mapping
#ace (ape) - maximum likelihood ancestral reconstruction
#
#
#
write.tree(ultra.tree, 'bok_ultra.tre')
ultra.tree <- ape::read.tree('bok_ultra.tre')



#select disc. variables
data_model_disc <- data %>% 
  select(base_da_aritenoide, apice_da_aritenoide, structure, emission_pattern, note_types) %>% 
  mutate(base_da_aritenoide  = case_when(data$base_da_aritenoide == "espesso" ~ 1,
                                         data$base_da_aritenoide == "espesso medio" ~ 2,
                                         data$base_da_aritenoide == "espesso com proeminencia" ~ 3)) %>% 
  mutate(apice_da_aritenoide = case_when(data$apice_da_aritenoide == "estreito" ~ 1,
                                         data$apice_da_aritenoide == "medio" ~ 2,
                                         data$apice_da_aritenoide == "ampla" ~ 3)) %>% 
  mutate(structure = ifelse(data$structure == 1, 2, 1)) %>% 
  mutate(emission_pattern = ifelse(data$emission_pattern == 1, 2, 1)) %>% 
  mutate(note_types = case_when(data$note_types == 2 ~ 3,
                                data$note_types == 1 ~ 2,
                                data$note_types == 0 ~ 1))
  
  
  



#build a treedata object and extract data with species names for each variable
data_model_disc <- treedata(ultra.tree, data_model_disc) %>% .$data




#model matrixes are defined in script analysis_discrete_evolmodels.R


##### base da aritenoide #####
  
  #stochastic character mapping
  #Q = "empirical" fits a Mk model to simulate all stochastic mapping (transition matrix); 
  #Q = "mcmc" computes nsim values of Q from the posterior MCMC bayesian simulations, then each simulated stochastic tree is generated for each Q matrix
  #Q = extracted from fitDiscrete from script "analysis_discrete_evolmodels.R"
  basearitenoide_bidir_stoc <- make.simmap(ultra.tree, data_model_disc[,'base_da_aritenoide'], model = 'ER', Q = as.Qmatrix(fit_ER_ba), nsim = 100)

  # extract q matrix
  q_basearitenoide <- as.Qmatrix(fit_ER_ba)
  plot(q_basearitenoide)
  
  # plot transition rates
  #tiff("rates_basearitenoide.tiff", units="in", width=8, height=8, res=300)
  plot(q_basearitenoide, rotate = 90, cex.rates=1, cex.traits=1, lwd = 3, spacer=0.20, text = TRUE)
  #dev.off()
  
  
  cols <- setNames(c("blue", "red", "light green"), levels(as.factor(data_model_disc[,"base_da_aritenoide"])))
  plot(summary(basearitenoide_bidir_stoc), offset = 0.5, ftype = "i", cex = 0.3, colors = cols)
  text(x=0.2, y=7, "Base da aritenoide", cex = 1)
  add.simmap.legend(leg = c("Espessamento sem proeminência", "Proeminência pouco desenvolvida", "Proeminência bem desenvolvida"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)


##### apice da aritenoide #####

  apicearitenoide_dir_stoc <- make.simmap(ultra.tree, data_model_disc[,'apice_da_aritenoide'], model = directional.model.3states.012, Q = as.Qmatrix(fit_directional_aa), nsim = 100)
  
  # extract q matrix
  q_apicearitenoide <- as.Qmatrix(fit_directional_aa)
  plot(q_apicearitenoide)
  
  # plot transition rates
  #tiff("rates_apicearitenoide.tiff", units="in", width=8, height=8, res=300)
  plot(q_apicearitenoide, rotate = 90, cex.rates=1, cex.traits=1, lwd = 3, spacer=0.20, text = TRUE)
  #dev.off()
  
  
  cols <- setNames(c("blue", "red", "light green"), levels(as.factor(data_model_disc[,"apice_da_aritenoide"])))
  plot(summary(apicearitenoide_dir_stoc), offset = 0.5, ftype = "i", cex = 0.3, colors = cols)
  text(x=0.2, y=7, "Ápice da aritenoide", cex = 1)
  add.simmap.legend(leg = c("Estreito", "Médio", "Amplo"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)
  
  


##### structure #####

  #run ancestral reconstruction with the best model chosen from fitDiscrete = directional model
  #ace (maximum likelihood): computes the maximum likelihood for each possible char state at nodes
  #struc_er <- ace(data_model_disc[,"structure"],ultra.tree,type="discrete",method="ML",model="ER")

  #get likelihoods of char states for each node
  #struc_er$lik.anc

  
  #stochastic character mapping
  #Q = "empirical" fits a Mk model to simulate all stochastic mapping (transition matrix); 
  #Q = "mcmc" computes nsim values of Q from the posterior MCMC bayesian simulations, then each simulated stochastic tree is generated for each Q matrix
  #Q = extracted from fitDiscrete from script "analysis_discrete_evolmodels.R"
  struc_dir_stoc <- make.simmap(ultra.tree, data_model_disc[,"structure"], model = directional.model, Q = as.Qmatrix(fit_directional_stru), nsim = 100)
  
  # extract q matrix
  q_struc <- as.Qmatrix(fit_directional_stru)
  plot(q_struc)

  # plot transition rates
  tiff("rates_struc.tiff", units="in", width=8, height=8, res=300)
  plot(q_struc, rotate = 90, cex.rates=1, cex.traits=1, lwd = 3, spacer=0.20, text = TRUE)
  dev.off()

#plot pies probabilities
  #from ace
  plotTree(ultra.tree, offset=0.5, ftype = 'i')
  nodelabels(pie = struc_ard$lik.anc, piecol = c("blue", "red"), cex = 0.4)

  
  #from make.simmap
  cols <- setNames(c("blue", "red"), levels(as.factor(data_model_disc[,"structure"])))
  plot(summary(struc_dir_stoc), offset = 0.5, ftype = "i", cex = 0.3, colors = cols)
  text(x=0.2, y=7, "Estrutura", cex = 1)
  add.simmap.legend(leg = c("Tonal", "Pulsado"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)


  #plot all stochastic mapped trees
  par(mfrow=c(10,10))
  sapply(struc_dir_stoc,plotSimmap,lwd=1,ftype="off")


  
##### emission pattern #####

  #run ancestral reconstruction with the best model chosen from fitDiscrete = directional model
  
  #ace
  #emission_er <- ace(data_model_disc[,"emission_pattern"],ultra.tree,type="discrete",method="ML",model="ER")
  #get likelihoods of char states for each node
  #emission_er$lik.anc
  
  
  #stochastic character mapping
  emission_dir_stoc <- make.simmap(ultra.tree, data_model_disc[,"emission_pattern"], model = directional.model, Q = as.Qmatrix(fit_directional_emis1), nsim = 100)
  
  # extract q matrix
  q_emis <- as.Qmatrix(fit_directional_emis1)
  plot(q_emis)
  
  # plot transition rates
  tiff("rates_emission_pattern.tiff", units="in", width=8, height=8, res=300)
  plot(q_emis, rotate = 90, cex.rates=1, cex.traits=1, lwd = 3, spacer=0.15, text = TRUE)
  dev.off()
  



  #plot pies probabilities
  #from ace
  plotTree(tree, offset=0.5)
  nodelabels(pie = emission_er$lik.anc, piecol = c("blue", "red"), cex = 0.4)
  
  
  #from make.simmap
  cols <- setNames(c("blue", "red"), levels(as.factor(data_model_disc[,"emission_pattern"])))
  plot(summary(emission_dir_stoc), offset = 0.5, ftype = "i", cex = 0.3, colors = cols)
  text(x=0.2, y=7, "Padrão de emissão de notas", cex = 1)
  add.simmap.legend(leg = c("Contínuo", "Série"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)

  
##### note_types ##### 

  #run ancestral reconstruction with the best model chosen from fitDiscrete = directional model
  
  #ace
  #note_er <- ace(data_model_disc[,"note_types"],ultra.tree,type="discrete",method="ML",model="ER")
  
  #get likelihoods of char states for each node
  #note_er$lik.anc
  
  
  #stochastic character mapping
  note_dir_stoc <- make.simmap(ultra.tree, data_model_disc[,"note_types"], model = directional.model1, Q = as.Qmatrix(fit_directional_note1), nsim = 100)
  
  # extract q matrix
  q_note <- as.Qmatrix(fit_directional_note1)
  plot(q_note)
  
  # plot transition rates
  tiff("rates_note_types.tiff", units="in", width=10, height=10, res=300)
  plot(q_note, rotate = 90, cex.traits='off', lwd = 5, spacer=0.35, text = FALSE)
  dev.off()
  
  

  
  
  #plot pies probabilities
  #from ace
  plotTree(tree, offset=0.5)
  nodelabels(pie = note_er$lik.anc, piecol = c("blue", "red", "green"), cex = 0.4)
  
  
  #from make.simmap
  cols <- setNames(c("blue", "red", "light green"), levels(as.factor(data_model_disc[,"note_types"])))
  plot(summary(note_dir_stoc), offset = 0.5, ftype = "i", cex = 0.3, colors = cols)
  text(x=0.2, y=7, "Tipos de notas", cex = 1)
  add.simmap.legend(leg = c("Um tipo", "Dois tipos (obrigatória)", "Dois tipos (facultativa)"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)
  
  
  
  
  
##### FINAL FIGURE #####

  #base da aritenoide
  tiff('basearitenoide_phylo.tiff', width = 22, height = 15, units = "cm", res = 300, compression = "lzw")
  
  cols <- setNames(c("blue", "red", "light green"), levels(as.factor(data_model_disc[,"base_da_aritenoide"])))
  plot(summary(basearitenoide_bidir_stoc), offset = 0.2, ftype = "i", cex = 0.5, colors = cols, fsize = 1.2)
  text(x=0.2, y=7, "Base da aritenoide", cex = 1)
  add.simmap.legend(leg = c("Sem proeminência (0)", "Proeminência pouco desenvolvida (1)", "Proeminência bem desenvolvida (2)"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)
  
  dev.off()
  
  #apice da aritenoide
  tiff('apicearitenoide_phylo.tiff', width = 22, height = 15, units = "cm", res = 300, compression = "lzw")
  
  
  cols <- setNames(c("blue", "red", "light green"), levels(as.factor(data_model_disc[,"apice_da_aritenoide"])))
  plot(summary(apicearitenoide_dir_stoc), offset = 0.2, ftype = "i", cex = 0.5, colors = cols, fsize = 1.2)
  text(x=0.2, y=7, "Ápice da aritenoide", cex = 1)
  add.simmap.legend(leg = c("Estreito (0)", "Médio (1)", "Amplo (2)"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)
  
  dev.off()
  
  
  #structure
  tiff('structure_phylo.tiff', width = 17, height = 15, units = "cm", res = 300, compression = "lzw")
  
  
  cols <- setNames(c("blue", "red"), levels(as.factor(data_model_disc[,"structure"])))
  plot(summary(struc_dir_stoc), offset = 0.2, ftype = "i", cex = 0.5, colors = cols, fsize = 1.2)
  text(x=0.2, y=7, "Estrutura", cex = 1)
  add.simmap.legend(leg = c("Tonal (0)", "Pulsado (1)"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)
  
  dev.off()
  
  #emission pattern
  tiff('emission_pattern_phylo.tiff', width = 17, height = 15, units = "cm", res = 300, compression = "lzw")
  
  
  cols <- setNames(c("blue", "red"), levels(as.factor(data_model_disc[,"emission_pattern"])))
  plot(summary(emission_dir_stoc), offset = 0.2, ftype = "i", cex = 0.5, colors = cols, fsize = 1.2)
  text(x=0.2, y=7, "Padrão de emissão", cex = 1)
  add.simmap.legend(leg = c("Contínuo (0)", "Série (1)"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)
  
  dev.off()
  
  #note types
  tiff('note_types_phylo.tiff', width = 17, height = 15, units = "cm", res = 300, compression = "lzw")
  
  cols <- setNames(c("blue", "red", "light green"), levels(as.factor(data_model_disc[,"note_types"])))
  plot(summary(note_dir_stoc), offset = 0.2, ftype = "i", cex = 0.5, colors = cols, fsize = 1.2)
  text(x=0.2, y=7, "Tipos de notas", cex = 1)
  add.simmap.legend(leg = c("Um tipo (0)", "Dois tipos (obrigatória) (1)", "Dois tipos (facultativa) (2)"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)
  
  dev.off()
  
