#Evaluating better evolutionary model fit through the fitContinuous function in geiger pckg.

#select cont. variables
data_model_cont <- data %>% select(call_duration, call_rate, call_freq_pk, call_freq_bw, mean_svl)
#z-scale data
data_model_cont <- scale(data_model_cont)
#build a treedata object and extract data with species names for each variable
data_model_cont <- treedata(tree, data_model_cont) %>% .$data


#create a vector of all evolutionary models possible
models <- c("BM","OU","EB","rate_trend","lambda",
            "kappa","delta","mean_trend","white")


##fitting the models for each variable and each evol. model
evol_model_cont <- apply(data_model_cont, 2, function(x) {
    sapply(models, function(i) {
      fit <- fitContinuous(ultra.tree, dat = x, model = i)
      
      data.frame(AIC = fit$opt$aic, 
                 AICc = fit$opt$aicc,
                 lnL = fit$opt$lnL)
  })
})



#likelihood ratio tests?






##########################Former approach##########################



#take each model and each variable in data_model and apply fitContinuous function
#extract AIC
evol_model_aic <- sapply(models, function(i) {
  apply(data_model, 2, function(x) {
    fit <- fitContinuous(tree, dat = x, model = i)
    fit$opt$aic
  })
})

##extract AICc
evol_model_aicc <- sapply(models, function(i) {
  apply(data_model, 2, function(x) {
    fit <- fitContinuous(tree, dat = x, model = i)
    fit$opt$aicc
  })
})

#extract lok likelihood
evol_model_lk <- sapply(models, function(i) {
  apply(data_model, 2, function(x) {
    fit <- fitContinuous(tree, dat = x, model = i)
    fit$opt$lnL
  })
})


#combine model fits into a single df
evol_model_aic <- t(evol_model_aic) %>% as.data.frame %>% 
  mutate(param = "AIC", .before = everything())

evol_model_aicc <- t(evol_model_aicc) %>% as.data.frame %>% 
  mutate(param = "AICc", .before = everything())

evol_model_lk <- t(evol_model_lk) %>% as.data.frame %>% 
  mutate(param = "lnLk", .before = everything())

results <- rbind(evol_model_aic, evol_model_aicc, evol_model_lk)
results


  
#compare different model performances visually for all traits
idx <- 0
apply(results, 2, function(trait, name) {
  idx <<- idx + 1;
  results %>% 
    ggplot(aes(x = rownames(results), y = trait)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(label = name[idx]) +
    xlab("Model") +
    ylab("") +
    facet_wrap(.~param, scale = "free")
}, colnames(results))







