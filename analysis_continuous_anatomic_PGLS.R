# inspect the volume and fibre length separately for exploratory analysis and see the variation inter-species
# check if muscles with similar volumes have different fibre lengths that can indicate longer fascicle length
# which can be inferred that those longer muscles will have more number of sarcomeres (larger contractile property)

# test for collinearity between muscle volumes and fibre lengths

# 
# 
# after fitting pgls model of svl vs pcsa, use residual values to fit anatomical traits vs. acoustic ones
# 
# use backwise step selection model to test best fit of anatomical traits vs categorical acoustic variables
# 
# explore pcsa functional space (PCSA plotted against fibre length)

# plot residuals from pgls for each species (e.g. larynx size vs body size)


  # take out the "chronos" from tree and import again only phylo to work with caper functions

  write.tree(ultra.tree, 'bok_ultra.tre')
  ultra.tree <- ape::read.tree('bok_ultra.tre')
  ggplot2::theme_set(ggplot2::theme_bw())
  
  #### ANATOMICAL TRAITS ####
  
  # wrangle, log transform and scale dataset to work with caper functions
  
  data_pgls <- data %>% 
    select(mean_svl, larynx_size:constrictor_post_pcsa) %>% 
    rename(svl = 'mean_svl') %>% 
    log() %>% 
    scale() 

  
  # make a species column to work with caper functions
  
  data_pgls <- cbind(species = rownames(data_pgls), 
                     as.data.frame(data_pgls)) 
  
  
  
  # create a comparative.data to compute pgls
  
  bok_pgls_data_an <- comparative.data(ultra.tree, data_pgls, species, vcv = TRUE)
  bok_pgls_data_an$data
  
  
  pgls_fit <- data.frame(formula = NA, lambda = NA, intercept = NA, coeff = NA, pvalue = NA, adj_r_squared = NA)
  residuals <- data.frame(row.names = rownames(bok_pgls_data_an$data))
  phy_residuals <- data.frame(row.names = rownames(bok_pgls_data_an$data))
  
  # PGLS size vs. all anatomical variables
  for (i in 1:(dim(bok_pgls_data_an$data)[2])) {
    
    if (colnames(bok_pgls_data_an$data)[i] == 'svl') {next}
    
    form <- formula(paste(colnames(bok_pgls_data_an$data)[i], ' ~ svl'))
    mod <- pgls(form, bok_pgls_data_an, lambda = 'ML')
    
    
    pgls_fit[i, 'formula'] <- format(form)
    pgls_fit[i, 'lambda'] <- summary(mod)$param[2]
    pgls_fit[i, 'intercept'] <- summary(mod)$coefficients[1,1]
    pgls_fit[i, 'coeff'] <- summary(mod)$coefficients[2,1]
    pgls_fit[i, 'pvalue'] <- summary(mod)$coefficients[2,4]
    pgls_fit[i, 'adj_r_squared'] <- summary(mod)$adj.r.squared
    
    residuals[, colnames(bok_pgls_data_an$data)[i]] <- mod$residuals
    phy_residuals[, colnames(bok_pgls_data_an$data)[i]] <- mod$phyres
  }
  

  pgls_fit <- pgls_fit %>% drop_na()
  pgls_fit
  #write_xlsx(pgls_fit, 'pgls_fit_anatomic.xlsx')
  
  ## RESIDUALS ##
  residuals
  phy_residuals
  
  

  
  
  ### fibre length vs muscle volumes ###

  
  # #pgls_fit for fibre lenght vs musc. vols
  # pgls_fit_fl <- data.frame(formula = NA, intercept = NA, coeff = NA, pvalue = NA, adj_r_squared = NA)
  # residuals_fl <- data.frame(row.names = rownames(bok_pgls_data_an$data))
  # 
  # formulas <- c('dilatator_fl ~ dilatator_vol', 'constrictor_ext_fl ~ constrictor_ext_vol',
  #               'constrictor_ant_fl ~ constrictor_ant_vol', 'constrictor_post_fl ~ constrictor_post_vol')
  #   
  # for(i in 1:length(formulas)){
  # 
  #   form <- formula(formulas[i])
  #   mod <- pgls(form, bok_pgls_data_an)
  #   
  #   
  #   pgls_fit_fl[i, 'formula'] <- format(form)
  #   pgls_fit_fl[i, 'intercept'] <- round(summary(mod)$coefficients[1,1],3)
  #   pgls_fit_fl[i, 'coeff'] <- round(summary(mod)$coefficients[2,1], 3)
  #   pgls_fit_fl[i, 'pvalue'] <- round(summary(mod)$coefficients[2,4], 5)
  #   pgls_fit_fl[i, 'adj_r_squared'] <- round(summary(mod)$adj.r.squared, 2)
  #   
  #   residuals_fl[, format(form)] <- mod$residuals
  # 
  # }
  # 
  # 
  #   pgls_fit_fl
  # 
  #   residuals_fl <- residuals_fl %>% rename(mdl = 1,
  #                                     mce = 2,
  #                                     mca = 3,
  #                                     mcp = 4)
  # 
  #   residuals_fl
  
  
  
    
    
  ## Test using phyl.resid from phytools ##
  ## 
  test_resid <- treedata(ultra.tree, data_pgls[,-1]) %>% .$data
  
  phytools_residuals <- phyl.resid(ultra.tree, x = test_resid[,'svl'], Y = test_resid[,'vc_area'])
  
  phytools_residuals$resid[,1]
  
  plotTree.barplot(ultra.tree, phytools_residuals$resid[,1])
  
  plot(phytools_residuals$resid[,1], phy_residuals$vc_area)
  ## phyl.resid matches the residuals from pgls and not phyresiduals
  ##
  
  
  
  
  
  
  #### PGLS PLOTS ####
  
  ## plot svl vs morphometric larynx traits
  
  # get line coeffs
  lines_ <- data.frame(intercept = pgls_fit[1:2, 'intercept'], 
                       slope = pgls_fit[1:2, 'coeff'], 
                       trait = factor(c('LL', 'VCA')))
  
  
  # make plot
  plot_pgls_morphometrics <- ggplot(bok_pgls_data_an$data) +
    geom_point(aes(svl, larynx_size), alpha = 0.7, size = 3) +
    geom_point(aes(svl, vc_area), shape = 0, alpha = 0.7, size = 3) +
    geom_abline(lines_, mapping = aes(slope = slope, intercept = intercept, linetype = trait)) +
    theme_bw() +
    theme(legend.position = c(.13,.85), #legend position are x and y offsets from bottom left
          legend.box.background = element_rect(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)
          ) +
    scale_linetype_discrete(name = "Trait", ) +
    labs(x = 'SVL (mm)', y = 'Larynx length and vocal cord area')
    

  # save fig
  tiff("plot_pgls_anatomic_morphometrics.tiff", units="cm", width=16, height=10, res=300)
  plot_pgls_morphometrics
  dev.off()  
  
  
  ## plot svl vs mdl architecture
  
  # melt data for ggplot
  dilatator_data <- bok_pgls_data_an$data %>% 
    select(svl, dilatator_vol, dilatator_fl, dilatator_pcsa) %>% 
    reshape2::melt(id='svl')
  
  
  # get line coeffs
  lines_dilatator = data.frame(intercept = pgls_fit[c(3,10,11), 'intercept'], 
                               slope = pgls_fit[c(3,10,11), 'coeff'], 
                               trait = factor(c('Volume', 'Comp. fibras', 'PCSA')))
  
  
  # make plot
  plot_pgls_dilatator <- ggplot(dilatator_data) +
    geom_point(aes(svl, value, shape = variable), alpha = 0.7, size = 3, show.legend = FALSE) +
    geom_abline(lines_dilatator, mapping = aes(slope = slope, intercept = intercept, linetype = trait)) +
    theme(
          legend.position = c(.15,.80), #legend position are x and y offsets from bottom left
          legend.box.background = element_rect(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)
          ) +
    scale_shape_manual("m. d. l.", values = c(16, 0, 2)) +
    scale_linetype_manual("m. d. l.", values = c("solid", "longdash", "dotted")) +
    labs(x = 'CRC (mm)', y = 'Volume, Comp. Fibras e PCSA')
    #scale_linetype_discrete(name = "m. d. l.", ) +
    #labs(linetype = "m. d. l.", shape = "m. d. l.") +
    
  
  
  # save fig
  tiff("plot_pgls_dilatator.tiff", units="cm", width=16, height=10, res=300)
  plot_pgls_dilatator
  dev.off()  
  
  
  
  # plot svl vs. mce architecture
  
  # melt data for ggplot
  mce_data <- bok_pgls_data_an$data %>% 
    select(svl, constrictor_ext_vol, constrictor_ext_fl, constrictor_ext_pcsa) %>% 
    reshape2::melt(id='svl')
  
  
  # get line coeffs
  lines_mce = data.frame(intercept = pgls_fit[c(4,8,12), 'intercept'], 
                               slope = pgls_fit[c(4,8,12), 'coeff'], 
                               trait = factor(c('Volume', 'Comp. fibras', 'PCSA')))
  
  # make plot
  plot_pgls_mce <- ggplot(mce_data) +
    geom_point(aes(svl, value, shape = variable), alpha = 0.7, size = 3, show.legend = FALSE) +
    geom_abline(lines_mce, mapping = aes(slope = slope, intercept = intercept, linetype = trait)) +
    theme(
      legend.position = c(.15,.80), #legend position are x and y offsets from bottom left
      legend.box.background = element_rect(),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    scale_shape_manual("m. c. l. e.", values = c(16, 0, 2)) +
    scale_linetype_manual("m. c. l. e.", values = c("solid", "longdash", "dotted")) +
    labs(x = 'CRC (mm)', y = 'Volume, Comp. Fibras e PCSA')
  
  
  
  # save fig
  tiff("plot_pgls_mce.tiff", units="cm", width=16, height=10, res=300)
  plot_pgls_mce
  dev.off()  
  
  
  
  # plot svl vs. mca architecture
  
  # melt data for ggplot
  mca_data <- bok_pgls_data_an$data %>% 
    select(svl, constrictor_ant_vol, constrictor_ant_fl, constrictor_ant_pcsa) %>% 
    reshape2::melt(id='svl')
  
  
  # get line coeffs
  lines_mca = data.frame(intercept = pgls_fit[c(5,7,13), 'intercept'], 
                         slope = pgls_fit[c(5,7,13), 'coeff'], 
                         trait = factor(c('Volume', 'Comp. fibras', 'PCSA')))
  
  
  # make plot
  plot_pgls_mca <- ggplot(mca_data) +
    geom_point(aes(svl, value, shape = variable), alpha = 0.7, size = 3, show.legend = FALSE) +
    geom_abline(lines_mca, mapping = aes(slope = slope, intercept = intercept, linetype = trait)) +
    theme(
      legend.position = c(.15,.80), #legend position are x and y offsets from bottom left
      legend.box.background = element_rect(),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    scale_shape_manual("m. c. l. a.", values = c(16, 0, 2)) +
    scale_linetype_manual("m. c. l. a.", values = c("solid", "longdash", "dotted")) +
    labs(x = 'CRC (mm)', y = 'Volume, Comp. Fibras e PCSA')
  
  
  
  # save fig
  tiff("plot_pgls_mca.tiff", units="cm", width=16, height=10, res=300)
  plot_pgls_mca
  dev.off()  
  
  
  
  # plot svl vs. mcp architecture
  
  # melt data for ggplot
  mcp_data <- bok_pgls_data_an$data %>% 
    select(svl, constrictor_post_vol, constrictor_post_fl, constrictor_post_pcsa) %>% 
    reshape2::melt(id='svl')
  
  
  # get line coeffs
  lines_mcp = data.frame(intercept = pgls_fit[c(6,9,14), 'intercept'], 
                         slope = pgls_fit[c(6,9,14), 'coeff'], 
                         trait = factor(c('Volume', 'Comp. fibras', 'PCSA')))
  
  
  # make plot
  plot_pgls_mcp <- ggplot(mcp_data) +
    geom_point(aes(svl, value, shape = variable), alpha = 0.7, size = 3, show.legend = FALSE) +
    geom_abline(lines_mcp, mapping = aes(slope = slope, intercept = intercept, linetype = trait)) +
    theme(
      legend.position = c(.15,.80), #legend position are x and y offsets from bottom left
      legend.box.background = element_rect(),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    scale_shape_manual("m. c. l. p.", values = c(16, 0, 2)) +
    scale_linetype_manual("m. c. l. p.", values = c("solid", "longdash", "dotted")) +
    labs(x = 'CRC (mm)', y = 'Volume, Comp. Fibras e PCSA')
  
  
  
  # save fig
  tiff("plot_pgls_mcp.tiff", units="cm", width=16, height=10, res=300)
  plot_pgls_mcp
  dev.off()  
  
  
  #grid arrange and save fig
  require(gridExtra)
  tiff("plot_pgls_muscles_svl.tiff", units="cm", width=30, height=20, res=300)
  grid.arrange(plot_pgls_dilatator, plot_pgls_mce, plot_pgls_mca, plot_pgls_mcp)
  dev.off()
  
  
  
  
