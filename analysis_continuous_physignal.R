#Anatomy
##Phylogenetic signal - calculate for the entire anatomical dataset?
##Phylogenetic signal - calculate for each continuous variable of larynx

#Acoustics
#Phylogenetic signal - calculate for the entire acoustic dataset?
#Phylogenetic signal - calculate for each continuous variable except pulse rates and temp

#select cont. variables
  data_physig_cont <- data %>% 
    select(call_duration, call_rate, call_freq_pk, call_freq_bw, mean_svl)

#z-scale data
  data_physig_cont <- scale(data_physig_cont)

#build a treedata object and extract data with species names for each variable
  data_physig_cont <- treedata(ultra.tree, data_physig_cont) %>% .$data

#build a standard error object
  se_physig_cont <- data %>% select(se_dur, se_rate, se_pkfreq, se_bwfreq, se_svl)


  
### Without measurement error
  
#use t() to transpose arrays to a matrix like structure
#use simplify2array to get arrays from results of apply
#apply() to apply phylosig function to each variable

#Pagel lambda
  lambda <- t(simplify2array(apply(data_physig_cont, #data obj
                         2,  #column axis
                         phylosig, #function and arguments
                         tree = ultra.tree, 
                         method = "lambda", 
                         test = T))) %>% 
    as.data.frame() %>% 
    select(lambda, P) %>% 
    rename("lambda_Pvalue" = "P")

#K Blomberg
  k_blomberg <- t(simplify2array(apply(data_physig_cont,
                         2,  
                         phylosig, 
                         tree = ultra.tree, 
                         method = "K", 
                         test = T))) %>% 
    as.data.frame() %>% 
    select(K, P) %>% 
    rename("K_Pvalue" = "P")


#compare both estimates of lambda and K
phylosignal_continuous <- bind_cols(lambda, k_blomberg)
phylosignal_continuous


### With measurement error

#solve standard errors NA














