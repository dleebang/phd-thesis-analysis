# get ultra.tree
  write.tree(ultra.tree, 'bok_ultra.tre')
  ultra.tree <- ape::read.tree('bok_ultra.tre')

  
# get relevant acoustic traits
  data_ac <- data %>% 
    select(!c(A_pulse_rate, B_pulse_rate, air_temp) & call_duration:call_freq_bw) %>% 
    log() %>% 
    scale() 
  
  data_ac <- cbind(species = rownames(data_ac), 
                   as.data.frame(data_ac), data[,c('structure', 'emission_pattern', 'note_types')]) 

  # data_ac['structure'] <- ifelse(data_ac['structure'] == 0, 'Tonal', 'Pulsado') %>% as.factor()
  # data_ac['emission_pattern'] <- ifelse(data_ac['emission_pattern'] == 0, 'Continuo', 'Serie') %>% as.factor()
  # data_ac['note_types'] <- case_when(data_ac['note_types'] == 0 ~ 'Um',
  #                                     data_ac['note_types'] == 1 ~ 'Dois_obrigatorio',
  #                                     data_ac['note_types'] == 2 ~ 'Dois_facultativo') %>% as.factor()


  
# get residuals from anatomic pgls analysis
  residuals_an <- residuals %>% select(larynx_size, vc_area, dilatator_pcsa:constrictor_post_pcsa)

  
# merge both datasets
  data_pgls <- cbind(as.matrix(residuals_an), data_ac, data[,c("base_da_aritenoide","apice_da_aritenoide")])
  
  
# make comparative data
  data_pgls <- comparative.data(ultra.tree, data_pgls, species, vcv = TRUE)
  
  data_pgls
  
  
# create dataframe to store results
  pgls_fit <- data.frame(formula = NA, lambda = NA, intercept = NA, coeff = NA, pvalue = NA, adj_r_squared = NA)


# run pgls for continuous traits (anatomic vs acoustic)
  for (i in colnames(data_pgls$data[c(1:6)])) {
    for (j in colnames(data_pgls$data[7:10])) {

    
    form <- formula(paste(j, ' ~ ', i)) # y (response) ~ x (predictor)
    mod <- pgls(form, data_pgls, lambda = 'ML')
    
    
    pgls_fit[nrow(pgls_fit) +1, 'formula'] <- format(form)
    pgls_fit[nrow(pgls_fit), 'lambda'] <- summary(mod)$param[2]
    pgls_fit[nrow(pgls_fit), 'intercept'] <- summary(mod)$coefficients[1,1]
    pgls_fit[nrow(pgls_fit), 'coeff'] <- summary(mod)$coefficients[2,1]
    pgls_fit[nrow(pgls_fit), 'pvalue'] <- summary(mod)$coefficients[2,4]
    pgls_fit[nrow(pgls_fit), 'adj_r_squared'] <- summary(mod)$adj.r.squared
    
   }
}
  
  pgls_fit <- pgls_fit %>% drop_na()
  pgls_fit
  apply(pgls_fit[,2:6], 2, function(x) round(x, 3))
  
  write_xlsx(pgls_fit, 'pgls_fit_anatomic_and_acoustic.xlsx')
  
  
  
  
  ##########
  # running pgls of discrete anatomic predictors with continuous acoustic tratis
  
  # merge both datasets
  data_pgls <- cbind(data_ac, data[,c("base_da_aritenoide","apice_da_aritenoide")])
  
  
  # make comparative data
  data_pgls <- comparative.data(ultra.tree, data_pgls, species, vcv = TRUE)
  
  data_pgls
  
  
  
  mod_disc_gls_base_aritenoide <- gls(call_freq_bw ~ base_da_aritenoide,
                                 data=data_pgls$data,correlation=
                                   corBrownian(1,ultra.tree),method="ML")
  
  summary(mod_disc_gls_base_aritenoide)
  anova(mod_disc_gls_base_aritenoide)
  
  
  
  mod_disc_gls_apice_aritenoide <- gls(call_freq_bw ~ apice_da_aritenoide,
                                      data=data_pgls$data,correlation=
                                        corBrownian(1,ultra.tree),method="ML")
  
  summary(mod_disc_gls_apice_aritenoide)
  anova(mod_disc_gls_apice_aritenoide)

  
  ##########
  
  
  
  
  
  
  ##########
  
  #running pgls of combined larynx traits to predict acoustic continuous traits
  
  #emission_pattern
  mod_disc_gls_em_pattern <- gls(emission_pattern ~ larynx_size + vc_area + 
                              dilatator_pcsa + 
                              constrictor_ext_pcsa + 
                              constrictor_ant_pcsa + 
                              constrictor_post_pcsa +
                                base_da_aritenoide +
                                apice_da_aritenoide,
                            data=data_pgls$data,correlation=
                              corBrownian(1,ultra.tree),method="ML")

  
  summary(mod_disc_gls_em_pattern)
  anova(mod_disc_gls_em_pattern)
  
  #stepAIC(mod_disc_gls_em_pattern, direction = 'backward', trace = 1)
  
  
  #note_types 
  mod_disc_gls_note_types <- gls(note_types ~ larynx_size + vc_area + 
                                   dilatator_pcsa + 
                                   constrictor_ext_pcsa + 
                                   constrictor_ant_pcsa + 
                                   constrictor_post_pcsa +
                                   base_da_aritenoide +
                                   apice_da_aritenoide,
                                 data=data_pgls$data,correlation=
                                   corBrownian(1,ultra.tree),method="ML")
  
  summary(mod_disc_gls_note_types)
  anova(mod_disc_gls_note_types)
  
  #stepAIC(mod_disc_gls_note_types, direction = 'backward', trace = 1)
  
  
 
  
  
  
  
  ### PLOTS
  ggplot2::theme_set(ggplot2::theme_bw())
  #significant correlations: pgls_fit 2, 5-8 (acv), 10, 24, 27, 29, 32
  pgls_fit
  data_pgls$data
  
  
  plot_pgls_2 <- ggplot(data_pgls$data, aes(larynx_size, call_rate)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[2,'coeff'], intercept = pgls_fit[2,'intercept']) +
    xlab('Comprimento da laringe (resíduos)') +
    ylab('T. emissão cantos') +
    theme(text = element_text(size = 15))  
  
  plot_pgls_5 <- ggplot(data_pgls$data, aes(call_duration, vc_area)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[5,'coeff'], intercept = pgls_fit[5,'intercept']) +
    xlab('Área da corda vocal (resíduos)') +
    ylab('Dur. cantos') +
    theme(text = element_text(size = 15))  
         
  plot_pgls_6 <- ggplot(data_pgls$data, aes(call_rate, vc_area)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[6,'coeff'], intercept = pgls_fit[6,'intercept']) +
    xlab('Área da corda vocal (resíduos)') +
    ylab('T. emissão cantos') +
    theme(text = element_text(size = 15))  
  
  plot_pgls_7 <- ggplot(data_pgls$data, aes(call_freq_pk, vc_area)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[7,'coeff'], intercept = pgls_fit[7,'intercept']) +
    xlab('Área da corda vocal (resíduos)') +
    ylab('Freq. dominante') +
    theme(text = element_text(size = 15))  
  
  plot_pgls_8 <- ggplot(data_pgls$data, aes(call_freq_bw, vc_area)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[8,'coeff'], intercept = pgls_fit[8,'intercept']) +
    xlab('Área da corda vocal (resíduos)') +
    ylab('B. Freq.') +
    theme(text = element_text(size = 15))  
  
  plot_pgls_10 <- ggplot(data_pgls$data, aes(call_rate, dilatator_pcsa)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[10,'coeff'], intercept = pgls_fit[10,'intercept']) +
    xlab('PCSA m. d. l.') +
    ylab('T. emissão cantos') +
    theme(text = element_text(size = 15))  
  
  plot_pgls_24 <- ggplot(data_pgls$data, aes(call_freq_bw, constrictor_post_pcsa)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[24,'coeff'], intercept = pgls_fit[24,'intercept']) +
    xlab('PCSA m. c. l. p.') +
    ylab('Banda de frequência') +
    theme(text = element_text(size = 15))  
  

  tiff("correlation_plots_anat_acoustic.tiff", units="in", width=16, height=10, res=300, compression = 'lzw')
  cowplot::plot_grid(plot_pgls_2,
                          plot_pgls_5,
                          plot_pgls_6,
                          plot_pgls_7,
                          plot_pgls_8,
                          plot_pgls_10,
                          plot_pgls_24,
                          labels="AUTO",
                     label_fontface = 'plain',
                     label_size = 22)
  dev.off()
  
  
  
  #plots for significant discrete predictors (base da aritenoide x freq. dominante, apice da aritenoide x dur canto, apice x banda de freq.)
  
  data_pgls$data[,'base_da_aritenoide'] <- ifelse(data_pgls$data[,'base_da_aritenoide'] == 'espesso medio', 'proeminência pouco desenvolvida', data_pgls$data[,'base_da_aritenoide'])
  data_pgls$data[,'base_da_aritenoide'] <- ifelse(data_pgls$data[,'base_da_aritenoide'] == 'espesso com proeminencia', 'proeminência bem desenvolvida', data_pgls$data[,'base_da_aritenoide'])
  
  
  tiff("boxplots_base_aritenoide.tiff", units="in", width=8, height=6, res=300, compression = 'lzw')
  
  ggplot(data_pgls$data, aes(x = base_da_aritenoide, y = call_freq_pk)) + 
    geom_boxplot() + 
    geom_point(alpha = 0.5, size = 3) +
    xlab('Base da aritenoide') +
    ylab('Frequência dominante')
    theme(text = element_text(size = 18))  
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  # larynx_traits <- c('larynx_size',
  #                    'larynx_size + vc_area',
  #                    'larynx_size + vc_area + dilatator_pcsa',
  #                    'larynx_size + vc_area + dilatator_pcsa + constrictor_ext_pcsa',
  #                    'larynx_size + vc_area + dilatator_pcsa + constrictor_ext_pcsa + constrictor_ant_pcsa',
  #                    'larynx_size + vc_area + dilatator_pcsa + constrictor_ext_pcsa + constrictor_ant_pcsa + constrictor_post_pcsa')
  # 
  # #call dur
  # pgls_fit_combined_dur <- list()
  # AIC_mod_dur <- data.frame(AIC = NA, AICc = NA, lik = NA)
  # 
  # for(i in 1:length(larynx_traits)) {
  #   formula <- formula(paste('call_duration ~ ', larynx_traits[i]))
  #   mod <- pgls(formula, data = data_pgls, lambda = 'ML')
  #   pgls_fit_combined_dur[[i]] <- summary(mod)
  #   names(pgls_fit_combined_dur)[i] <- larynx_traits[i]
  #   AIC_mod_dur[i, 'AIC'] <- mod$aic
  #   AIC_mod_dur[i, 'AICc'] <- mod$aicc
  #   AIC_mod_dur[i, 'lik'] <- mod$model$log.lik
  #   rownames(AIC_mod_dur)[i] <- larynx_traits[i]
  # }
  # 
  # pgls_fit_combined_dur
  # setNames(aic.w(AIC_mod_dur$AIC), larynx_traits)
  # 
  # 
  # #call rate
  # pgls_fit_combined_rate <- list()
  # AIC_mod_rate <- data.frame(AIC = NA, AICc = NA, lik = NA)
  # 
  # for(i in 1:length(larynx_traits)) {
  #   formula <- formula(paste('call_rate ~ ', larynx_traits[i]))
  #   mod <- pgls(formula, data = data_pgls, lambda = 'ML')
  #   pgls_fit_combined_rate[[i]] <- summary(mod)
  #   names(pgls_fit_combined_rate)[i] <- larynx_traits[i]
  #   AIC_mod_rate[i, 'AIC'] <- mod$aic
  #   AIC_mod_rate[i, 'AICc'] <- mod$aicc
  #   AIC_mod_rate[i, 'lik'] <- mod$model$log.lik
  #   rownames(AIC_mod_rate)[i] <- larynx_traits[i]
  # }
  # 
  # pgls_fit_combined_rate
  # setNames(aic.w(AIC_mod_rate$AIC), larynx_traits)
  # 
  # 
  # #pk freq
  # pgls_fit_combined_pkfreq <- list()
  # AIC_mod_pkfreq <- data.frame(AIC = NA, AICc = NA, lik = NA)
  # 
  # for(i in 1:length(larynx_traits)) {
  #   formula <- formula(paste('call_freq_pk ~ ', larynx_traits[i]))
  #   mod <- pgls(formula, data = data_pgls, lambda = 'ML')
  #   pgls_fit_combined_pkfreq[[i]] <- summary(mod)
  #   names(pgls_fit_combined_pkfreq)[i] <- larynx_traits[i]
  #   AIC_mod_pkfreq[i, 'AIC'] <- mod$aic
  #   AIC_mod_pkfreq[i, 'AICc'] <- mod$aicc
  #   AIC_mod_pkfreq[i, 'lik'] <- mod$model$log.lik
  #   rownames(AIC_mod_pkfreq)[i] <- larynx_traits[i]
  # }
  # 
  # pgls_fit_combined_pkfreq
  # setNames(aic.w(AIC_mod_pkfreq$AIC), larynx_traits)
  # 
  # 
  # #freq bw
  # pgls_fit_combined_freqbw <- list()
  # AIC_mod_freqbw <- data.frame(AIC = NA, AICc = NA, lik = NA)
  # 
  # for(i in 1:length(larynx_traits)) {
  #   formula <- formula(paste('call_duration ~ ', larynx_traits[i]))
  #   mod <- pgls(formula, data = data_pgls, lambda = 'ML')
  #   pgls_fit_combined_freqbw[[i]] <- summary(mod)
  #   names(pgls_fit_combined_freqbw)[i] <- larynx_traits[i]
  #   AIC_mod_freqbw[i, 'AIC'] <- mod$aic
  #   AIC_mod_freqbw[i, 'AICc'] <- mod$aicc
  #   AIC_mod_freqbw[i, 'lik'] <- mod$model$log.lik
  #   rownames(AIC_mod_freqbw)[i] <- larynx_traits[i]
  # }
  # 
  # pgls_fit_combined_freqbw
  # setNames(aic.w(AIC_mod_freqbw$AIC), larynx_traits)
  # 
  # 
  # 
  # 
  # #running pgls of larynx traits to predict discrete traits
  # #structure
  # mod_disc_gls_struc <- gls(structure ~ larynx_size + vc_area + 
  #                             dilatator_pcsa + 
  #                             constrictor_ext_pcsa + 
  #                             constrictor_ant_pcsa + 
  #                             constrictor_post_pcsa +
  #                             base_da_aritenoide +
  #                             apice_da_aritenoide,
  #                           data=data_pgls$data,correlation=
  #                             corBrownian(1,ultra.tree),method="ML")
  # 
  # summary(mod_disc_gls_struc)
  # anova(mod_disc_gls_struc)
  # 
  # #stepAIC(mod_disc_gls_struc, direction = 'backward', trace = 1)
  