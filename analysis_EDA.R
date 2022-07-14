#Average and sd
avg_sd <- bok_traits %>% group_by(species) %>% 
          summarize(Average = mean(call_duration, na.rm = TRUE), 
          Std_deviation = sd(call_duration, na.rm = TRUE))

###boxplots
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 13)) #set center title for all ggplots


## call duration
p_calldur <- bok_traits %>% filter(!is.na(call_duration)) %>%
   mutate(species = reorder(species, call_duration, FUN = median)) %>%
   ggplot(aes(call_duration, species)) 

p_calldur + geom_boxplot() + scale_x_continuous(trans = "log10") +
   geom_vline(aes(xintercept=median(call_duration))) +
   ggtitle("Duração do canto (s)") + 
   labs(x = "", y  = "")


## call rate
p_rate <- bok_traits %>% filter(!is.na(call_rate)) %>%
   mutate(species = reorder(species, call_rate, FUN = median)) %>%
   ggplot(aes(call_rate, species))

p_rate + geom_boxplot() + scale_x_continuous(trans = "log10") +
  geom_vline(aes(xintercept=median(call_rate))) +
  ggtitle("Taxa de emissão de cantos (cantos/min)") + 
  labs(x = "", y  = "")
  

## call pk frequency
p_pkfreq <- bok_traits %>% filter(!is.na(call_freq_pk)) %>%
  mutate(species = reorder(species, call_freq_pk, FUN = median)) %>%
  ggplot(aes(call_freq_pk, species))

p_pkfreq + geom_boxplot() + scale_x_continuous(trans = "log10") +
  geom_vline(aes(xintercept=median(call_freq_pk))) +
  ggtitle("Frequência dominante do canto (Hz)") + 
  labs(x = "", y  = "")


## call freq bwdith
p_bwfreq <- bok_traits %>% filter(!is.na(call_freq_bw)) %>%
  mutate(species = reorder(species, call_freq_bw, FUN = median)) %>%
  ggplot(aes(call_freq_bw, species))

p_bwfreq + geom_boxplot() + scale_x_continuous(trans = "log10") +
  geom_vline(aes(xintercept=median(call_freq_bw))) +
  ggtitle("Banda de frequência do canto (Hz)") +
  labs(x = "", y  = "")


## svl
p_svl <- bok_svl %>% filter(!is.na(svl)) %>%
  mutate(species = reorder(species, svl, FUN = median)) %>%
  ggplot(aes(svl, species))

p_svl + geom_boxplot() + scale_x_continuous(trans = "log10") +
  geom_vline(aes(xintercept=median(svl))) +
  ggtitle("Comprimento rostro-cloacal (mm)") +
  labs(x = "", y  = "")



###histograms
data_hist <- data %>% select(mean_svl,
                             call_duration, 
                             call_rate, 
                             call_freq_pk, 
                             call_freq_bw)


#with Rbase
par(mfrow = c(2,3))
par(mar=c(4, 4.5, 4, 4.5))
for (i in 1:ncol(data_hist)) {
  trait_name <- parse_character(colnames(data_hist)[i])
  hist(data_hist[,i],
       prob = T,
       main = trait_name,
       xlab = "")
  lines(density(data_hist[,i]))
}


#with ggplot
data_hist %>% gather(key = trait) %>% 
  ggplot(aes(value)) +
  geom_bar(alpha = 0.8) +
  scale_x_binned() +
  theme_bw() +
  xlab("") +
  facet_wrap(.~trait, scale = "free")




### scatterplot
data %>% filter(!is.na(call_freq_pk) & !is.na(mean_svl)) %>%
  ggplot(aes(mean_svl, call_freq_pk)) +
  geom_point(alpha = 0.6)

cor(data$mean_svl, data$call_freq_pk)

