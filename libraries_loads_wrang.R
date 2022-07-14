#Load libraries
libraries <- c("caper", 
               "geiger", 
               "ape", 
               "phytools",
               "phangorn",
               "phylobase",
               "readxl", 
               "castor", 
               "tidyverse", 
               "nlme", 
               "gtools")


if(sum(as.numeric(!libraries %in% installed.packages())) != 0){
  installer <- libraries[!libraries %in% installed.packages()]
  for(i in 1:length(installer)) {
    install.packages(installer, dependencies = T)
    break()}
  sapply(libraries, require, character = T) 
} else {
  sapply(libraries, require, character = T) 
}

#change sci notation to standard notation
options(scipen=0)

#Load tree
tree <- read.tree(file = "Boktree_newick.txt")

#Ultrametricize tree
ultra.tree <- chronos(tree, lambda = 0)
is.ultrametric(ultra.tree)


#Load traits dataset
bok_traits <- read_excel("Bok_calltraits.xlsx", col_types = c("text",
                                                              "text",
                                                              "text",
                                                              "text",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric",
                                                              "numeric"))

#Load SVL dataset
bok_svl <- read_excel("Bok_svl.xlsx")

#inspect dataset and tree
head(bok_traits)
str(bok_traits)
str(tree)

#Wrangling:

  #calculate mean and se for trait dataset
  data <- bok_traits %>% 
    group_by(species, clade) %>% 
    mutate(se_dur = sd(call_duration, na.rm=T)/sqrt(n()),
           se_rate = sd(call_rate, na.rm=T)/sqrt(n()),
           se_pkfreq = sd(call_freq_pk, na.rm=T)/sqrt(n()),
           se_bwfreq = sd(call_freq_bw, na.rm=T)/sqrt(n())) %>% 
    summarize_if(is.numeric, mean, na.rm=T) %>% 
    ungroup()
  
  #missing value imputation for B. vulcaniae with mean values of its clade
  #call rate
  data$call_rate[is.na(data$call_rate)] <- data %>% 
    filter(clade == "circumdata") %>% 
    summarize(cr_mean = mean(call_rate, na.rm=T)) %>% 
    .$cr_mean

  
  #call_freq_bw
  data$call_freq_bw[is.na(data$call_freq_bw)] <- data %>% 
    filter(clade == "circumdata") %>% 
    summarize(bw_mean = mean(call_freq_bw, na.rm=T)) %>% 
    .$bw_mean
  
  
  #calculate mean for svl dataset
  svl <- bok_svl %>%
    group_by(species) %>%
    mutate(se_svl = sd(svl, na.rm=T)/sqrt(n())) %>% 
    summarize(mean_svl = mean(svl, na.rm = T), se_svl = mean(se_svl, na.rm = T))
  
  
  #merge data and svl
  data <- data %>% 
    inner_join(svl)
  
  
  #make species names to be rownames
  data <- data %>% 
    remove_rownames %>% 
    column_to_rownames(var="species")
  
  
  #keep only epithet names in tree
  ultra.tree$tip.label <- str_replace(ultra.tree$tip.label, ".*Bok_", "") %>% 
    str_replace("\\_.*", "")
  
  
  #Check if tip labels match species name in the dataset
  check <- name.check(ultra.tree, data)
  
  
  #Drop tip with no data
  ultra.tree <- drop.tip(ultra.tree, check$tree_not_data)

  
  #Sort data rows in the same order as they are in tree$tip.label
  data <- data[order(match(rownames(data), ultra.tree$tip.label)),]
  identical(rownames(data), ultra.tree$tip.label)
  
plotTree(ultra.tree)
