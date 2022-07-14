#compute ancestral state reconstructions based on best fit model (ARD, SYM, ER)
#compute ancestral state reconstructions based on all discrete evol models and compare likelihoods?
#
#
#make.simmap (phytools) - stochastic character mapping
#ace (ape) - maximum likelihood ancestral reconstruction
#
#

#select disc. variables
data_model_disc <- data %>% 
  select(structure, emission_pattern, note_types)


#build a treedata object and extract data with species names for each variable
data_model_disc <- treedata(ultra.tree, data_model_disc) %>% .$data


##### structure #####

#run ancestral reconstruction with the best model chosen from fitDiscrete
  #ace (maximum likelihood): computes the maximum likelihood for each possible char state at nodes
  struc_ard <- ace(data_model_disc[,"structure"],ultra.tree,type="discrete",method="ML",model="ARD")
  
  #get likelihoods of char states for each node
  struc_ard$lik.anc

  
  #stochastic character mapping
  #Q = "empirical" uses most likely value for Q (transition matrix); 
  #Q = "mcmc" computes nsim values of Q from the posterior MCMC bayesian simulations, then each simulated stochastic tree is generated for each Q matrix
  struc_ard_bayes <- make.simmap(ultra.tree, data_model_disc[,"structure"], model = "ARD", Q = "empirical", nsim = 100)

  describe.simmap(struc_ard_bayes)

  
#plot pies probabilities
  #from ace
  plotTree(ultra.tree, offset=0.5)
  nodelabels(pie = struc_ard$lik.anc, piecol = c("blue", "red"), cex = 0.4)

  
  #from make.simmap
  cols <- setNames(c("blue", "red"), levels(as.factor(data_model_disc[,"structure"])))
  plot(summary(struc_ard_bayes), offset = 0.5, ftype = "i", cex = 0.3, colors = cols)
  text(x=0.2, y=7, "Structure", cex = 1)
  add.simmap.legend(leg = c("Tonal", "Pulsed"),
                    colors = cols,
                    prompt = FALSE,
                    x = 0.15, y = 5,
                    fsize = 0.85)


  #plot all stochastic mapped trees
  par(mfrow=c(10,10))
  sapply(struc_ard_bayes,plotSimmap,lwd=1,ftype="off")


  
##### emission pattern #####
emission_er <- ace(data_model_disc[,"emission_pattern"],tree,type="discrete",method="ML",model="ER")
emission_er_bayes <- make.simmap(tree, data_model_disc[,"emission_pattern"], model = "ER", Q = "mcmc", nsim = 100)


#get likelihoods of char states for each node
emission_er$lik.anc


#plot pies probabilities
#from ace
plotTree(tree, offset=0.5)
nodelabels(pie = emission_er$lik.anc, piecol = c("blue", "red"), cex = 0.4)


#from make.simmap
plot(summary(emission_er_bayes), offset = 0.5, ftype = "i", cex = 0.3)


##### note_types ##### 
note_er <- ace(data_model_disc[,"note_types"],tree,type="discrete",method="ML",model="ER")
note_er_bayes <- make.simmap(tree, data_model_disc[,"note_types"], model = "ER", Q = "mcmc", nsim = 100)


#get likelihoods of char states for each node
note_er$lik.anc


#plot pies probabilities
#from ace
plotTree(tree, offset=0.5)
nodelabels(pie = note_er$lik.anc, piecol = c("blue", "red", "green"), cex = 0.4)


#from make.simmap
plot(summary(note_er_bayes), offset = 0.5, ftype = "i", cex = 0.3)
