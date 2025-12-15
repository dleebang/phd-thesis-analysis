#Phylogenetic signal - calculate for the entire discrete dataset?
#Phylogenetic signal - calculate for each discrete variable

#phylo.sig.disc function (Maddison & Slatkin p value)
source("Func.phylo.signal.disc.txt")

#delta statistics from Borges et al. 2019
source("delta_statistics.R")

# get ultra.tree
write.tree(ultra.tree, 'bok_ultra.tre')
ultra.tree <- ape::read.tree('bok_ultra.tre')


#### ANATOMY
#select disc. anatomic variables
data_physig <- data %>% 
  select(base_da_aritenoide, apice_da_aritenoide) %>% 
  mutate(base_da_aritenoide  = case_when(data$base_da_aritenoide == "espesso" ~ 0,
                   data$base_da_aritenoide == "espesso medio" ~ 1,
                   data$base_da_aritenoide == "espesso com proeminencia" ~ 2)) %>% 
  mutate(apice_da_aritenoide = case_when(data$apice_da_aritenoide == "estreito" ~ 0,
                                         data$apice_da_aritenoide == "medio" ~ 1,
                                         data$apice_da_aritenoide == "ampla" ~ 2))


#build a treedata object and extract data with species names for each variable
data_physig <- treedata(ultra.tree, data_physig) %>% .$data


#### Maddison & Slatkin
#run function for all discrete traits
results <- apply(data_physig, 2, phylo.signal.disc, ultra.tree, rep = 500)


#create a matrix to hold results
phylosignal_madslat <- matrix(ncol = 1, nrow = length(results), 
                              dimnames = list(trait = names(results), "P.value"))


#extract p values from results
for (i in names(results)) {
  phylosignal_madslat[i,] <- results[[i]]$.Randomization.Results["P.value",]
}

phylosignal_madslat

#### Delta statistic
#calculate delta and p-values
random_delta <- rep(NA,500)
phylosignal_delta <- apply(data_physig, 2, function(trait) {
  delta_trait <- delta(trait,ultra.tree,0.1,0.5,10000,10,100)
  for (i in 1:length(random_delta)){
    rtrait <- sample(trait)
    random_delta[i] <- delta(rtrait,ultra.tree,0.1,0.5,10000,10,100)
  }
  data.frame(delta_value = delta_trait, 
             p_value = sum(random_delta>delta_trait)/length(random_delta))
})

phylosignal_delta





#### ACOUSTIC
#select disc. acoustic variables
data_physig <- data %>% 
  select(structure, emission_pattern, note_types) %>% 
  scale()


#build a treedata object and extract data with species names for each variable
data_physig <- treedata(ultra.tree, data_physig) %>% .$data

#character states coding:
#structure
#   0 = tonal
#   1 = pulsed
#
#emission_pattern
#   0 = continuous
#   1 = well defined series
#
#note_types
#   0 = one
#   1 = two obligatory
#   2 = two facultative



#### Maddison & Slatkin
#run function for all discrete traits
results <- apply(data_physig, 2, phylo.signal.disc, ultra.tree, rep = 500)


#create a matrix to hold results
phylosignal_madslat <- matrix(ncol = 1, nrow = length(results), 
                          dimnames = list(trait = names(results), "P.value"))


#extract p values from results
for (i in names(results)) {
  phylosignal_madslat[i,] <- results[[i]]$.Randomization.Results["P.value",]
}

phylosignal_madslat





#### Delta statistic
#calculate delta and p-values
random_delta <- rep(NA,500)
phylosignal_delta <- apply(data_physig, 2, function(trait) {
  delta_trait <- delta(trait,ultra.tree,0.1,0.5,10000,10,100)
  for (i in 1:length(random_delta)){
    rtrait <- sample(trait)
    random_delta[i] <- delta(rtrait,ultra.tree,0.1,0.5,10000,10,100)
  }
  data.frame(delta_value = delta_trait, 
             p_value = sum(random_delta>delta_trait)/length(random_delta))
})

phylosignal_delta



##calculate delta for each trait 500 times to build confidence interval
##DONT RUN
delta_values <- t(replicate(500, apply(data_physig, 2, delta, ultra.tree,
                     lambda0 = 0.1, 
                     se = 0.5,
                     sim = 10000,
                     thin = 10,
                     burn = 100)))

head(delta_values)

CI_delta <- apply(delta_values, 2, function(x) 
  {X <- mean(x)
  SE <- sd(x) / sqrt(length(x))
  data.frame(lower = X - (qnorm(.975) * SE), 
             upper = X + (qnorm(.975) * SE))
  })

CI_delta
##






