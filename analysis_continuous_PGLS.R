# Correlated evolution between traits  

# PGLS assuming no phylogenetic covariance (=OLS)
# PGLS assuming Brownian motion (covariance represents proportional branch lenghts)
# PGLS with incorporated phylogenetic signal (lambda = 1 or 0) on VCV matrix (this is important due to sample size! Kamilar & Cooper 2013)
# PGLS calculate estimates of both regression parameters of two traits and lambda (phy. signal) on residuals of a regression

# PGLS for multimodel inference according to AICc values (backwise step model selection)
# https://stats.stackexchange.com/questions/35353/why-applying-model-selection-using-aic-gives-me-non-significant-p-values-for-the
# https://stats.stackexchange.com/questions/212674/model-selection-in-pgls
# report both p-value approach for individual variables and AICc values of backwise step model selection



  # take out the "chronos" from tree and import again only phylo to work with caper functions
  
  write.tree(ultra.tree, 'bok_ultra.tre')
  ultra.tree <- ape::read.tree('bok_ultra.tre')


  #### ACOUSTIC TRAITS ####
  
  # wrangle, log transform and scale dataset to work with caper functions
  
  data_pgls <- data %>% 
    select(mean_svl, !c(A_pulse_rate, B_pulse_rate, air_temp) & call_duration:call_freq_bw) %>% 
    rename(svl = 'mean_svl') %>% 
    log() %>% 
    scale() 


  # make a species column to work with caper functions
  
  data_pgls <- cbind(species = rownames(data_pgls), 
                     as.data.frame(data_pgls)) 


  # create a comparative.data to compute pgls
  
  bok_pgls_data_ac <- comparative.data(ultra.tree, data_pgls, species, vcv = TRUE)
  bok_pgls_data_ac$data

  
  pgls_fit <- data.frame(formula = NA, lambda = NA, intercept = NA, coeff = NA, pvalue = NA, adj_r_squared = NA)

  
  # PGLS size vs. all acoustic variables

  for (i in 2:(dim(bok_pgls_data_ac$data)[2])) {
    
      if (colnames(bok_pgls_data_ac$data)[i] == 'svl') {next}
      
      form <- formula(paste(colnames(bok_pgls_data_ac$data)[i], ' ~ svl'))
      mod <- pgls(form, bok_pgls_data_ac, lambda = 'ML')
      
  
      pgls_fit[i, 'formula'] <- format(form)
      pgls_fit[i, 'lambda'] <- summary(mod)$param[2]
      pgls_fit[i, 'intercept'] <- summary(mod)$coefficients[1,1]
      pgls_fit[i, 'coeff'] <- summary(mod)$coefficients[2,1]
      pgls_fit[i, 'pvalue'] <- summary(mod)$coefficients[2,4]
      pgls_fit[i, 'adj_r_squared'] <- summary(mod)$adj.r.squared
    }

  
  #PGLS call traits
  
  ## call duration vs. all variables
  
  for (i in 1:(dim(bok_pgls_data_ac$data[-(1:2)])[2])) {
    
    if (colnames(bok_pgls_data_ac$data[-(1:2)])[i] == 'call_duration') {next}
    
    form <- formula(paste('call_duration ~ ', colnames(bok_pgls_data_ac$data[-(1:2)])[i]))
    mod <- pgls(form, bok_pgls_data_ac, lambda = 'ML')
    
    
    pgls_fit[nrow(pgls_fit) + 1, 'formula'] <- format(form)
    pgls_fit[nrow(pgls_fit), 'lambda'] <- summary(mod)$param[2]
    pgls_fit[nrow(pgls_fit), 'intercept'] <- summary(mod)$coefficients[1,1]
    pgls_fit[nrow(pgls_fit), 'coeff'] <- summary(mod)$coefficients[2,1]
    pgls_fit[nrow(pgls_fit), 'pvalue'] <- summary(mod)$coefficients[2,4]
    pgls_fit[nrow(pgls_fit), 'adj_r_squared'] <- summary(mod)$adj.r.squared
  }

  
  ## call rate vs. pk freq. and bw freq.
  
  for (i in 1:(dim(bok_pgls_data_ac$data[-(1:3)])[2])) { #-(1:2) indexing removes species, svl and call dur from dataframe
    
    if (colnames(bok_pgls_data_ac$data[-(1:3)])[i] == 'call_rate') {next}
    
    form <- formula(paste('call_rate ~ ', colnames(bok_pgls_data_ac$data[-(1:3)])[i]))
    mod <- pgls(form, bok_pgls_data_ac, lambda = 'ML')
    
    
    pgls_fit[nrow(pgls_fit) + 1, 'formula'] <- format(form)
    pgls_fit[nrow(pgls_fit), 'lambda'] <- summary(mod)$param[2]
    pgls_fit[nrow(pgls_fit), 'intercept'] <- summary(mod)$coefficients[1,1]
    pgls_fit[nrow(pgls_fit), 'coeff'] <- summary(mod)$coefficients[2,1]
    pgls_fit[nrow(pgls_fit), 'pvalue'] <- summary(mod)$coefficients[2,4]
    pgls_fit[nrow(pgls_fit), 'adj_r_squared'] <- summary(mod)$adj.r.squared
  }
  
  
  ## pk freq. vs bw freq.
  
  mod <- pgls(call_freq_pk ~ call_freq_bw, bok_pgls_data_ac, lambda = 'ML')
  
  pgls_fit[nrow(pgls_fit) + 1, 'formula'] <- 'call_freq_pk ~ call_freq_bw'
  pgls_fit[nrow(pgls_fit), 'lambda'] <- summary(mod)$param[2]
  pgls_fit[nrow(pgls_fit), 'intercept'] <- summary(mod)$coefficients[1,1]
  pgls_fit[nrow(pgls_fit), 'coeff'] <- summary(mod)$coefficients[2,1]
  pgls_fit[nrow(pgls_fit), 'pvalue'] <- summary(mod)$coefficients[2,4]
  pgls_fit[nrow(pgls_fit), 'adj_r_squared'] <- summary(mod)$adj.r.squared
  
  
  ## drop NAs in pgls_fit
  
  pgls_fit <- pgls_fit %>% drop_na()
  pgls_fit
  write_xlsx(pgls_fit, 'pgls_fit_acoustic.xlsx')
  
  
  
  ### pgls 4, 5, 6, 8
  ggplot2::theme_set(ggplot2::theme_bw())

  plot_pgls_4 <- ggplot(bok_pgls_data_ac$data, aes(svl, call_freq_pk)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[3,'coeff'], intercept = pgls_fit[3,'intercept']) +
    xlab('CRC') +
    ylab('Frequência dominante') +
    theme(text = element_text(size = 20))      
  
  
  plot_pgls_5 <- ggplot(bok_pgls_data_ac$data, aes(call_rate, call_duration)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[5,'coeff'], intercept = pgls_fit[5,'intercept']) +
    xlab('Taxa de emissão de cantos') +
    ylab('Duração dos cantos') +
    theme(text = element_text(size = 20))  
  
  
  plot_pgls_6 <- ggplot(bok_pgls_data_ac$data, aes(call_freq_pk, call_duration)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[6,'coeff'], intercept = pgls_fit[6,'intercept']) +
    xlab('Frequência dominante') +
    ylab('Duração dos cantos') +
    theme(text = element_text(size = 20))  
  
  
  plot_pgls_8 <- ggplot(bok_pgls_data_ac$data, aes(call_freq_pk, call_rate)) +
    geom_point(size = 3) +
    geom_abline(slope = pgls_fit[8,'coeff'], intercept = pgls_fit[8,'intercept']) +
    xlab('Frequência dominante') +
    ylab('Taxa de emissão de cantos') +
    theme(text = element_text(size = 20))  
  
  
  tiff("correlation_plots_acoustic.tiff", units="in", width=16, height=10, res=300, compression = 'lzw')
  cowplot::plot_grid(plot_pgls_4,
                     plot_pgls_5,
                     plot_pgls_6,
                     plot_pgls_8,
                     labels="AUTO",
                     label_fontface = 'plain',
                     label_size = 22)
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #### CORRELOGRAM OF ALL TRAITS ####
  
  
  # plot correlogram of mean svl vs. call variables
  library(GGally)
  ggplot2::theme_set(ggplot2::theme_bw())
  
  bok_pgls_data_ac$data <- bok_pgls_data_ac$data[,-1]
  
  gg <- ggpairs(bok_pgls_data_ac$data, upper = list(continuous = 'blank'), 
                columnLabels = c('CRC',                                 
                                 'Duração',
                                 'Cantos/min',
                                 'F. dominante',
                                 'F. banda'))
  
  
  # set up GGally instance to customize plot to include pgls lines and coeffs
  
  m_ <- GGally::print_if_interactive
  

  # create column combinations and select those with variables names used to compute pgls
  
  cols <- expand.grid(names(bok_pgls_data_ac$data)[1:5], names(bok_pgls_data_ac$data)[1:4], stringsAsFactors=F)    
  cols <- cols[c(2:5, 8:10, 14:15, 20),] 
  
  
  # create a list to store plots to be later included in ggpairs
  
  plots_low <- list() #for lower diagonal
  plots_upper <- list() #for upper diagonal
  
  # get plot coordinates
  
  cols_right <- c(2:5, 3:5, 4:5, 5)
  cols_left <- rep(c(1, 2, 3, 4), times=c(4,3,2,1))
  
  
  # loop through cols object to create a customized ggplot with pgls lines, intercept and coeffs
  
  for (i in 1:nrow(cols)){
     p <- ggplot(bok_pgls_data_ac$data, aes_string(x = cols[i,2], y = cols[i,1]))+
          geom_point() +
            geom_abline(intercept = pgls_fit[i,]$intercept, slope = pgls_fit[i,]$coeff)
     plots_low[[i]] <- p
          
  }
  
  # check
  
  plots_low
  
  
  # loop through generated plots to insert customized plots within ggpairs lower diagonal
  
  for (i in seq_along(plots_low)){
    gg <- putPlot(gg, plots_low[[i]], cols_right[i], cols_left[i])
  }


  # now a loop to make upper diagonal plots with coeff and p values
  
  for (i in 1:nrow(pgls_fit)){
    p <- m_(ggally_text(paste('Coef: ', pgls_fit[i,]$coeff, '\n',
                              'p: ', pgls_fit[i,]$pvalue, '\n',
                              'R2 adj: ', pgls_fit[i,]$adj_r_squared), mapping = ggplot2::aes(color = I('black')))) + 
      theme(panel.grid.major = element_blank(), panel.border = element_blank())
    plots_upper[[i]] <- p
    
  }
 
  # check
  
  plots_upper

  
  # loop through generated plots to insert customized plots within ggpairs upper diagonal
  
  for (i in seq_along(plots_upper)){
    gg <- putPlot(gg, plots_upper[[i]], cols_left[i], cols_right[i]) #remember to invert coordinates
  }
  
  
  # make final adjustments
  gg[2,3] <- m_(ggally_text(paste('Coef: ', pgls_fit[5,]$coeff, '\n',
                                  'p: ', '<0.00001', '\n',
                                  'R2 adj: ', pgls_fit[5,]$adj_r_squared), mapping = ggplot2::aes(color = I('black')))) + 
    theme(panel.grid.major = element_blank(), panel.border = element_blank())
  
  
  # plot ggpairs
  
  gg
  
  # save img
  tiff("correlogram-acoustic.tiff", units="in", width=9, height=7, res=300, compression = 'lzw')
  gg
  dev.off()
  
  
  # code adapted from: https://stackoverflow.com/questions/30858337/how-to-customize-lines-in-ggpairs-ggally
  
  ###################################
  
  #### One per one plot approach ####
  #  m[1,2] <- m_(ggally_text(paste('Coef: ', pgls_fit[1,]$coeff, '\n',
  #'p: ', pgls_fit[1,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #  theme(panel.grid.major = element_blank())
  # 
  # m[1,3] <- m_(ggally_text(paste('Coef: ', pgls_fit[2,]$coeff, '\n',
  #                                'p: ', pgls_fit[2,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[1,4] <- m_(ggally_text(paste('Coef: ', pgls_fit[3,]$coeff, '\n',
  #                                'p: ', pgls_fit[3,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[1,5] <- m_(ggally_text(paste('Coef: ', pgls_fit[4,]$coeff, '\n',
  #                                'p: ', pgls_fit[4,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[2,3] <- m_(ggally_text(paste('Coef: ', pgls_fit[5,]$coeff, '\n',
  #                                'p: ', '<0.00001'), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[2,4] <- m_(ggally_text(paste('Coef: ', pgls_fit[6,]$coeff, '\n',
  #                                'p: ', pgls_fit[6,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[2,5] <- m_(ggally_text(paste('Coef: ', pgls_fit[7,]$coeff, '\n',
  #                                'p: ', pgls_fit[7,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[3,4] <- m_(ggally_text(paste('Coef: ', pgls_fit[8,]$coeff, '\n',
  #                                'p: ', pgls_fit[8,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[3,5] <- m_(ggally_text(paste('Coef: ', pgls_fit[9,]$coeff, '\n',
  #                                'p: ', pgls_fit[9,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # m[4,5] <- m_(ggally_text(paste('Coef: ', pgls_fit[10,]$coeff, '\n',
  #                                'p: ', pgls_fit[10,]$pvalue), mapping = ggplot2::aes(color = I('black')))) + 
  #   theme(panel.grid.major = element_blank())
  # 
  # 
  # m[2,1] <- m[2,1] + geom_abline(intercept = pgls_fit[1,]$intercept, slope = pgls_fit[1,]$coeff) 
  # m[3,1] <- m[3,1] + geom_abline(intercept = pgls_fit[2,]$intercept, slope = pgls_fit[2,]$coeff) 
  # m[4,1] <- m[4,1] + geom_abline(intercept = pgls_fit[3,]$intercept, slope = pgls_fit[3,]$coeff) 
  # m[5,1] <- m[5,1] + geom_abline(intercept = pgls_fit[4,]$intercept, slope = pgls_fit[4,]$coeff) 
  # m[3,2] <- m[3,2] + geom_abline(intercept = pgls_fit[5,]$intercept, slope = pgls_fit[5,]$coeff) 
  # m[4,2] <- m[4,2] + geom_abline(intercept = pgls_fit[6,]$intercept, slope = pgls_fit[6,]$coeff) 
  # m[5,2] <- m[5,2] + geom_abline(intercept = pgls_fit[7,]$intercept, slope = pgls_fit[7,]$coeff) 
  # m[4,3] <- m[4,3] + geom_abline(intercept = pgls_fit[8,]$intercept, slope = pgls_fit[8,]$coeff)
  # m[5,3] <- m[5,3] + geom_abline(intercept = pgls_fit[9,]$intercept, slope = pgls_fit[9,]$coeff)
  # m[5,4] <- m[5,4] + geom_abline(intercept = pgls_fit[10,]$intercept, slope = pgls_fit[10,]$coeff)
  # 
  # m  
