########
# Meta-analysis
########
## Packages
# Load libraries
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans, ggdist, rotl, ape)
# devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)

## Data
# bring in data and arrange for analysis
data <- read.csv("./Final.Analysis.Data/Meta_extracted_data.csv")
data$obs <- 1:nrow(data)
  
# Calculate effect sizes
data <- escalc(measure = "ROM", m1=trt_mean, sd1=trt_sd, n1=trt_n, m2=c_mean, sd2=c_sd , n2=c_n, data = data, append = TRUE)
escalc(measure = "ROM", m1=10, sd1=2, n1=10, m2=5, sd2=2 , n2=10)

# Models
# Intercept model, all random
model <- rma.mv(yi ~ 1, V = vi, random = list(~1|study_ID, ~1| species, ~1|obs), data = data)
saveRDS(model, "./Final.Models/Meta_analysis_models/Meta_mods/meta.model.rds")
# M1
model1 <- rma.mv(yi ~ trait,  random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model1, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m1.rds")
# M2
model2 <- rma.mv(yi ~ 0 + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model2, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m2.rds")
# M3
model3 <- rma.mv(yi ~ 0 + trait + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model3, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m3.rds")
# M4
model4 <- rma.mv(yi ~ 1 + age + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model4, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m4.rds")
# M5
model5 <- rma.mv(yi ~ 1 + scale(latitude, scale=FALSE) + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model5, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m5.rds")
# M6
model6 <- rma.mv(yi ~ 0 + class + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model6, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m6.rds")
# M7
model7 <- rma.mv(yi ~ 1 + species, random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model7, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m7.rds")
# M8
model8 <- rma.mv(yi ~ 0 + zone + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model8, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m8.rds")
# M9: intercept-based model
model9 <- rma.mv(yi ~ 1 + trait+ scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model9, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m9.rds")
# M10
model10 <- rma.mv(yi ~ 0 + habitat + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
saveRDS(model10, "./Final.Models/Meta_analysis_models/Meta_mods/meta.m10.rds")

# Funnel plots
# M2 - publication bias  
funnel( yaxis="seinv", model2, legend=FALSE, digit = 2, las = 1)
# log ratio of means
funnel( yaxis="seinv", model, digit = 2, las = 1)

# calculating reported I^2 measures from MCMCglmm model 
# devtools::install_github("daniel1noble/metaAidR")
library(metaAidR)
I2(model, v = data$vi, phylo = FALSE, obs = "obs")
I2(model3, v = data$vi, phylo = FALSE, obs = "obs")
I2(model4, v = data$vi, phylo = FALSE, obs = "obs")
I2(model5, v = data$vi, phylo = FALSE, obs = "obs")
I2(model6, v = data$vi, phylo = FALSE, obs = "obs")
I2(model8, v = data$vi, phylo = FALSE, obs = "obs")
I2(model7, v = data$vi, phylo = FALSE, obs = "obs")
  

### P. Pottier - meta-analysis with acclimation response ratio #### 
data <- data %>%  mutate(ARR=((trt_mean-c_mean)/(trt_temp-c_temp)),
                           Var_ARR = ((1/(trt_temp-c_temp))^2*(c_sd^2/c_n+trt_sd^2/trt_n)))
# acclimation response ratio figure
  ggplot(data) + geom_vline(xintercept=0, alpha=0.5, linetype=2)+
    stat_slab(aes(x = ARR, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.5,0.8, 0.95), labels = scales::percent_format()))), side = "bottom", scale = 0.5,
              show.legend = F, col = "darkcyan") + 
    stat_dots(aes(x = ARR, col=trait), alpha = 0.8, quantiles = 54,dotsize = 0.8, shape = 16)

# phylo tree - arranging data for 
  data<-as.data.frame(data)
  data$genus_species<-paste(data$genus, data$species)
  taxa<-tnrs_match_names(names=unique(data$genus_species), context="Animals") # match species names with the OTL taxonomy
  phylo_tree<-tol_induced_subtree(ott_ids=taxa$ott_id, label_format="name") # generate tree
  phylo_tree<-compute.brlen(phylo_tree, method="Grafen", power=1) # generate branch lengths
  # plot
  plot(phylo_tree)
# arranging figure
  phylo_tree$tip.label<-gsub("_", " ", phylo_tree$tip.label) # remove underscores in species names
  phylo_matrix<-vcv(phylo_tree, cor=T) # generate phylogenetic vcv matrix
  
  
##### Run models #####
# Model ALL: Both Tpref and CTmax analysed together, with all random effects
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
  # check phylogeny  
  orchaRd::i2_ml(model_all_rand) #  phylogeny explains virtually no variance, so we can remove it. 
  saveRDS(model_all_rand, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.phylo.rds")
  
# Intercept model: without phylogeny
  int_model<-rma.mv(ARR~1, 
                            V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                            method="REML",
                            test="t",
                            dfs="contain",
                            random=list(~1|study_ID, 
                                        ~1|obs),# Consider related_temps / related_individuals
                            data=data)
  saveRDS(int_model, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.no.phylo.rds")
  # model checks Intercept model
  summary(int_model)     
  i2_ml(int_model) # Lots of heterogeneity


# Trait model: Trait differences
  model_trait<- rma.mv(ARR~trait-1, 
                       V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|obs),# Consider related_temps / related_individuals
                       data=data)
  saveRDS(model_trait, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.trait.rds")
  summary(model_trait)
  orchard_plot(model_trait, group = "trait", mod="trait", xlab="ARR", data = data)
  
  
# Model species: differences between species
  model_species<- rma.mv(ARR~genus_species-1, 
                       V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|obs),# Consider related_temps / related_individuals
                       data=data)
  saveRDS(model_species, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.spp.rds")
  summary(model_species)
  orchard_plot(model_species, group = "genus_species", mod="genus_species", xlab="ARR", data = data)

  
## Life stage - AGE
  model_age<- rma.mv(ARR~age-1, 
                         V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                         method="REML",
                         test="t",
                         dfs="contain",
                     random=list(~1|study_ID, 
                                 ~1|obs),# Consider related_temps / related_individuals
                     data=data)
  saveRDS(model_age, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.age.rds")
  summary(model_age)
  orchard_plot(model_age, group = "age", mod="age", xlab="ARR", data = data)

  
## Latitude 
  model_zone<- rma.mv(ARR~zone-1, 
                       V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",random=list(~1|study_ID, 
                                                 ~1|obs),# Consider related_temps / related_individuals
                      data=data)
  saveRDS(model_zone, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.zone.rds")
  summary(model_zone)
  orchard_plot(model_zone, group = "zone", mod="zone", xlab="ARR", data = data)
  

## taxon
  model_taxon<- rma.mv(ARR~class-1, 
                     V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                     method="REML",
                     test="t",
                     dfs="contain",
                     random=list(~1|study_ID, 
                                 ~1|obs),# Consider related_temps / related_individuals
                     data=data)
  saveRDS(model_taxon, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.taxon.rds")
  summary(model_taxon)
  orchard_plot(model_taxon, group = "class", mod="class", xlab="ARR", data = data)
  
  