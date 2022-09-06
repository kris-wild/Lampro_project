
# Load libraries
  #install orchard plots
install.packages("R.rsp")
install.packages("devtools")
  devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)
  pacman::p_load(readxl, tidyverse, metafor, orchaRd, patchwork)

 # Data
  data <- read.csv("./Meta-analyses/extracted_data.csv")
  data$obs <- 1:nrow(data)
  
# Calculate effect sizes
  data <- escalc(measure = "ROM", m1=trt_mean, sd1=trt_sd, n1=trt_n, m2=c_mean, sd2=c_sd , n2=c_n, data = data, append = TRUE)
  
  escalc(measure = "ROM", m1=10, sd1=2, n1=10, m2=5, sd2=2 , n2=10)
  escalc(measure = "ROM", m1=5, sd1=2, n1=10, m2=10, sd2=2 , n2=10)
# Models
  model <- rma.mv(yi ~ 1, V = vi, random = list(~1|study_ID, ~1| species, ~1|obs), data = data)
  
  model1 <- rma.mv(yi ~ trait,  random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  
  model2 <- rma.mv(yi ~ 0 + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)

<<<<<<< HEAD
# Funnel plot for publication bias
  funnel(model2)
  

=======
  model3 <- rma.mv(yi ~ 0 + trait + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  
  model4 <- rma.mv(yi ~ 1 + age + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  
  model5 <- rma.mv(yi ~ 1 + scale(latitude, scale=FALSE) + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  
  model6 <- rma.mv(yi ~ 0 + class + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  
  model7 <- rma.mv(yi ~ 1 + species, random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  
  model8 <- rma.mv(yi ~ 0 + zone + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  
  model10 <- rma.mv(yi ~ 0 + habitat + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  #intercept-based model
  model9 <- rma.mv(yi ~ 1 + trait+ scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
  


  funnel( yaxis="seinv", model2, legend=FALSE, digit = 2, las = 1)
  
  # Funnel plot for publication bias
  funnel( yaxis="seinv", model, digit = 2, las = 1)

  library(devtools)
  install_github("daniel1noble/metaAidR")
  
  I2(model, v = data$vi, phylo = FALSE, obs = "obs")
  
  I2(model3, v = data$vi, phylo = FALSE, obs = "obs")
  
  I2(model4, v = data$vi, phylo = FALSE, obs = "obs")
  
  I2(model5, v = data$vi, phylo = FALSE, obs = "obs")
  
  I2(model6, v = data$vi, phylo = FALSE, obs = "obs")
  
  I2(model8, v = data$vi, phylo = FALSE, obs = "obs")
  
  I2(model7, v = data$vi, phylo = FALSE, obs = "obs")
  
>>>>>>> d045c38af99a226f2f62fa10f1889f150f9635d3

  
  

  ### P. Pottier - meta-analysis with acclimation response ratio #### 
  
  data <- data %>%  mutate(ARR=((trt_mean-c_mean)/(trt_temp-c_temp)),
                           Var_ARR = ((1/(trt_temp-c_temp))^2*(c_sd^2/c_n+trt_sd^2/trt_n)))
  
library(ggdist)
  ggplot(data) + geom_vline(xintercept=0, alpha=0.5, linetype=2)+
    stat_slab(aes(x = ARR, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.5,0.8, 0.95), labels = scales::percent_format()))), side = "bottom", scale = 0.5,
              show.legend = F, col = "darkcyan") + 
    stat_dots(aes(x = ARR, col=trait), alpha = 0.8, quantiles = 54,dotsize = 0.8, shape = 16)
  
  data<-as.data.frame(data)
  
  library(rotl)
  library(ape)
  data$genus_species<-paste(data$genus, data$species)
  taxa<-tnrs_match_names(names=unique(data$genus_species), context="Animals") # match species names with the OTL taxonomy
  phylo_tree<-tol_induced_subtree(ott_ids=taxa$ott_id, label_format="name") # generate tree
  phylo_tree<-compute.brlen(phylo_tree, method="Grafen", power=1) # generate branch lengths
  plot(phylo_tree)
  
  phylo_tree$tip.label<-gsub("_", " ", phylo_tree$tip.label) # remove underscores in species names
  phylo_matrix<-vcv(phylo_tree, cor=T) # generate phylogenetic vcv matrix
  
  
  ##### Run models #####
  
  # Both Tpref and CTmax analysed together, with all random effects
  
  model_all_rand<-rma.mv(ARR~1, 
                         V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                         method="REML",
                         test="t",
                         dfs="contain",
                         random=list(~1|study_ID, 
                                     ~1|genus_species,
                                     ~1|obs),# Consider related_temps / related_individuals
                         R=list(genus_species=phylo_matrix),
                         data=data)
  i2_ml(model_all_rand) #  phylogeny explains virtually no variance, so we can remove it. 
  
  
  # Intercept model without phylogeny
  int_model<-rma.mv(ARR~1, 
                            V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                            method="REML",
                            test="t",
                            dfs="contain",
                            random=list(~1|study_ID, 
                                        ~1|obs),# Consider related_temps / related_individuals
                            data=data)
  summary(int_model)     
  i2_ml(int_model) # Lots of heterogeneity
  
  
  orchard_plot(int_model, mod="Int", xlab="ARR")
  
  
  
  # Differences between traits
  
  model_trait<- rma.mv(ARR~trait-1, 
                       V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|es_ID),# Consider related_temps / related_individuals
                       data=data)
  summary(model_trait)
  orchard_plot(model_trait, mod="trait", xlab="ARR")
  
  
  # Differences between species
  
  model_species<- rma.mv(ARR~genus_species-1, 
                       V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|es_ID),# Consider related_temps / related_individuals
                       data=data)
  summary(model_species)

  # Other moderators 
  
  ## Life stage 
  
  model_stage<- rma.mv(ARR~age-1, 
                         V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                         method="REML",
                         test="t",
                         dfs="contain",
                         random=list(~1|study_ID, 
                                     ~1|es_ID),# Consider related_temps / related_individuals
                         data=data)
  summary(model_stage)
  orchard_plot(model_age, mod="age", xlab="ARR")

  ## Latitude 
  
  model_zone<- rma.mv(ARR~zone-1, 
                       V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|es_ID),# Consider related_temps / related_individuals
                       data=data)
  summary(model_zone)
  orchard_plot(model_zone, mod="latitude", xlab="ARR")
  
  
  ## taxon
  
  model_taxon<- rma.mv(ARR~class-1, 
                     V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                     method="REML",
                     test="t",
                     dfs="contain",
                     random=list(~1|study_ID, 
                                 ~1|es_ID),# Consider related_temps / related_individuals
                     data=data)
  summary(model_taxon)
  orchard_plot(model_taxon, mod="class", xlab="ARR")
  
  