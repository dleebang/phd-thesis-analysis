#Phylogenetic signal - calculate for the entire discrete dataset?
#Phylogenetic signal - calculate for each discrete variable

#select disc. variables
data_physig <- data %>% 
  select(structure, emission_pattern, note_types)

#z-scale data
data_physig <- scale(data_physig)

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


#phylo.sig.disc function 
source("Func.phylo.signal.disc.txt")

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



#delta statistics from Borges et al. 2019
source("delta_statistics.R")


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






