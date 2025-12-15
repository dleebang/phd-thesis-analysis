#fitPagel method of phytools to fit Pagel's (1994) test of correlated evolution of binary traits
#to test dependent model vs independent model of evolution

write.tree(ultra.tree, 'bok_ultra.tre')
ultra.tree <- ape::read.tree('bok_ultra.tre')



#select disc. variables and binarize data
data_model_disc <- data %>% 
  select(base_da_aritenoide, apice_da_aritenoide, structure, emission_pattern, note_types) %>% 
  mutate(base_da_aritenoide  = case_when(data$base_da_aritenoide == 'espesso medio' ~ 'proeminencia', 
                                         data$base_da_aritenoide == 'espesso com proeminencia' ~ 'proeminencia',
                                         data$base_da_aritenoide == 'espesso' ~ 'espesso')) %>%
  mutate(apice_da_aritenoide = ifelse(data$apice_da_aritenoide == 'medio', 'estreito', data$apice_da_aritenoide)) %>% 
  mutate(structure = ifelse(data$structure == 1, 'tonal', 'pulsado')) %>% 
  mutate(emission_pattern = ifelse(data$emission_pattern == 1, 'serie', 'continuo')) %>% 
  mutate(note_types = case_when(data$note_types == 2 ~ 'dois tipos',
                                data$note_types == 1 ~ 'dois tipos',
                                data$note_types == 0 ~ 'um tipo'))


#convert cols to factor to plot data matrix
data_model_disc_plot <- data_model_disc %>% 
  mutate_if(is.character, as.factor)



#build a treedata
data_model_disc <- treedata(ultra.tree, data_model_disc) %>% .$data

## ANATOMIC VS ANATOMIC
fit_anatomic <- fitPagel(ultra.tree, data_model_disc[,'base_da_aritenoide'], 
                         data_model_disc[,'apice_da_aritenoide'], 
                         dep.var = 'x'
                        )


## ACOUSTICS VS ACOUSTICS
fit_structure_emispat <- fitPagel(ultra.tree, data_model_disc[,'structure'], 
                                  data_model_disc[,'emission_pattern'], 
                                  dep.var = 'x'
                                  )


fit_structure_notetype <- fitPagel(ultra.tree, data_model_disc[,'structure'], 
                data_model_disc[,'note_types'], 
                dep.var = 'x'
                )


fit_emispat_notetype <- fitPagel(ultra.tree, data_model_disc[,'emission_pattern'], 
                                   data_model_disc[,'note_types'], 
                                   dep.var = 'x'
                                  )


## ANATOMIC VS ACOUSTICS


fitPagel(ultra.tree, data_model_disc[,'base_da_aritenoide'], 
         data_model_disc[,'structure'], 
         dep.var = 'x'
)

fitPagel(ultra.tree, data_model_disc[,'base_da_aritenoide'], 
         data_model_disc[,'emission_pattern'], 
         dep.var = 'x'
)

fitPagel(ultra.tree, data_model_disc[,'base_da_aritenoide'], 
         data_model_disc[,'note_types'], 
         dep.var = 'x'
)


fitPagel(ultra.tree, data_model_disc[,'apice_da_aritenoide'], 
         data_model_disc[,'structure'], 
         dep.var = 'x'
)

fitPagel(ultra.tree, data_model_disc[,'apice_da_aritenoide'], 
         data_model_disc[,'emission_pattern'], 
         dep.var = 'x'
)

fitPagel(ultra.tree, data_model_disc[,'apice_da_aritenoide'], 
         data_model_disc[,'note_types'], 
         dep.var = 'x'
)









# make plots for correlated anatomic traits
  tiff("corr_discrete_anatomic.tiff", units="in", width=7, height=10, res=300, compression = 'lzw')
  
  object <- plotTree.datamatrix(ultra.tree, data_model_disc_plot[,1:2], header=F, palettes=c("YlOrRd","PuBuGn"), fsize = 1.2)

  leg<-legend(x="bottomleft",
              names(object$colors$base_da_aritenoide),cex=1.2,
              pch=22,pt.bg=object$colors$base_da_aritenoide,pt.cex=1.5,
              bty="n",title="base da aritenoide")
  
  leg<-legend(x=leg$rect$left, y=leg$rect$top+leg$rect$h,
              names(object$colors$apice_da_aritenoide),
              cex=1.2,pch=22,pt.bg=object$colors$apice_da_aritenoide,
              pt.cex=1.5,bty="n",title="apice da aritenoide")
  
  
  dev.off()
  
  
  
  
  tiff("rates_corr_models_anatomic.tiff", units="in", width=7, height=10, res=300, compression = 'lzw')
  plot(fit_anatomic, cex.sub=1.2, cex.traits=1, cex.rates = 1.1, lwd.by.rate=T)
  dev.off()




# make plots for correlated acoustic traits
  tiff("corr_discrete_notetypes_emispat.tiff", units="in", width=7, height=10, res=300, compression = 'lzw')
  
  object <- plotTree.datamatrix(ultra.tree, data_model_disc_plot[,4:5], header=F, palettes=c("YlOrRd","PuBuGn"), fsize = 1.2)
  
  leg<-legend(x="bottomleft", inset = 0.05,
              names(object$colors$note_types),cex=1.2,
              pch=22,pt.bg=object$colors$note_types,pt.cex=1.5,
              bty="n",title="tipos de notas")
  
  leg<-legend(x=leg$rect$left, y=leg$rect$top+leg$rect$h,
              names(object$colors$emission_pattern),
              cex=1.2,pch=22,pt.bg=object$colors$emission_pattern,
              pt.cex=1.5,bty="n",title="padrão de emissão")
  
  
  dev.off()
  
  
  
  
  tiff("rates_corr_models_acoustic.tiff", units="in", width=7, height=10, res=300, compression = 'lzw')
  plot(fit_emispat_notetype, lwd.by.rate = T, cex.sub=1.2, cex.traits=1, cex.rates = 1.1)
  dev.off()
  



















# http://blog.phytools.org/2014/12/r-function-for-pagels-1994-correlation.html
# http://blog.phytools.org/search?q=correlated+evolution+discrete+traits
# http://blog.phytools.org/2016/05/some-cool-additional-features-for.html
