########
# Meta-analysis
########
## Packages
# Load libraries
pacman::p_load(devtools, tidyverse, dplyr, metafor, patchwork, R.rsp, emmeans, ggdist, rotl, ape, orchaRd, patchwork, latex2exp, cowplot)

## data
# bring in data and arrange for analysis
data <- read.csv("./Final.Analysis.Data/Final_Meta_results_extracted_data.csv")
data$obs <- 1:nrow(data)

# Calculate effect sizes and check 
data <- escalc(measure = "ROM", m1=trt_mean, sd1=trt_sd, n1=trt_n, m2=c_mean, sd2=c_sd , n2=c_n, data = data, append = TRUE)
escalc(measure = "ROM", m1=10, sd1=2, n1=10, m2=5, sd2=2 , n2=10)
model.bias <- rma.mv(yi ~ 0 + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
# Funnel plot
funnel( yaxis="seinv", model.bias, legend=FALSE, digit = 2, las = 1)

### ### ### 
### Meta-analysis with acclimation response ratio #### 
### ### ### 
data <- data %>%  mutate(ARR=((trt_mean-c_mean)/(trt_temp-c_temp)),
                         Var_ARR = ((1/(trt_temp-c_temp))^2*(c_sd^2/c_n+trt_sd^2/trt_n)))
# acclimation response ratio figure
ggplot(data) + geom_vline(xintercept=0, alpha=0.5, linetype=2)+
  stat_slab(aes(x = ARR, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(0.5,0.8, 0.95), labels = scales::percent_format()))), side = "bottom", scale = 0.5,
            show.legend = F, col = "darkcyan") + 
  stat_dots(aes(x = ARR, col=trait), alpha = 0.8, quantiles = 54,dotsize = 0.8, shape = 16)

### ### ###  
# phylo tree 
### ### ### 
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

### ### ###  
##### All models & checks #####
### ### ### 

# 1) Model ALL: Both Tpref and CTmax analysed together, with all random effects
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

# Drop phylogeny and just check if species is important. Basically the same so keep species in.
model_all_rand_spp<-rma.mv(ARR~1, 
                       V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|genus_species,
                                   ~1|obs),# Consider related_temps / related_individuals
                       data=data)

# 2) Intercept model: without phylogeny
int_model<-rma.mv(ARR~1, 
                  V=Var_ARR,
                  method="REML",
                  test="t",
                  dfs="contain",
                  random=list(~1|study_ID, 
                              ~1|genus_species,
                              ~1|obs),# Consider related_temps / related_individuals
                  data=data)
saveRDS(int_model, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.no.phylo.rds")
# model checks Intercept model
summary(int_model)     
i2_ml(int_model) # Lots of heterogeneity
predict(int_model) # Prediction intervals

# 3) Trait model: Trait differences
model_trait<- rma.mv(ARR~trait-1, 
                     V=Var_ARR,                     method="REML",
                     test="t",
                     dfs="contain",
                     random=list(~1|study_ID,
                                 ~1|genus_species,
                                 ~1|obs),# Consider related_temps / related_individuals
                     data=data)
saveRDS(model_trait, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.trait.rds")
summary(model_trait)
orchard_plot(model_trait, group = "study_ID", mod="trait", xlab="ARR", data = data)
i2_ml(model_trait)

# 4) Model species: differences between species
model_species<- rma.mv(ARR~genus_species-1, 
                       V=Var_ARR,
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID,
                                   ~1|genus_species,
                                   ~1|obs),# Consider related_temps / related_individuals
                       data=data)
saveRDS(model_species, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.spp.rds")
summary(model_species)
orchard_plot(model_species, group = "study_ID", mod="genus_species", xlab="ARR", data = data)
i2_ml(model_species)

## 5) Life stage - AGE
model_age<- rma.mv(ARR~age-1, 
                   V=Var_ARR,
                   method="REML",
                   test="t",
                   dfs="contain",
                   random=list(~1|study_ID, 
                               ~1|genus_species,
                               ~1|obs),# Consider related_temps / related_individuals
                   data=data)
saveRDS(model_age, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.age.rds")
summary(model_age)
Lifestage_plot <- orchard_plot(model_age, group = "study_ID", mod="age", xlab="ARR", data = data, legend.pos = "none") + 
  scale_fill_manual(values = c("deepskyblue2","firebrick", "darkgoldenrod1"))+ 
  scale_color_manual(values = c("deepskyblue2","firebrick", "darkgoldenrod1"))+
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle = 45))
i2_ml(model_age)


# 6) Latitude 
model_zone<- rma.mv(ARR~zone-1, 
                    V=Var_ARR,
                    method="REML",
                    test="t",
                    dfs="contain",random=list(~1|study_ID, 
                                              ~1|genus_species,
                                              ~1|obs),# Consider related_temps / related_individuals
                    data=data)
saveRDS(model_zone, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.zone.rds")
summary(model_zone)
LatLong_plot <-orchard_plot(model_zone, group = "study_ID", mod="zone", xlab="ARR", data = data, legend.pos = "none") +
  scale_fill_manual(values = c("deepskyblue2","firebrick"))+ 
  scale_color_manual(values = c("deepskyblue2","firebrick"))+
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text.y = element_text(angle = 45))
i2_ml(model_zone)


## 7) taxon
model_taxon<- rma.mv(ARR~class-1, 
                     V=Var_ARR,# Consider using the VCV matrix of ARR with correlated errors
                     method="REML",
                     test="t",
                     dfs="contain",
                     random=list(~1|study_ID, 
                                 ~1|genus_species,
                                 ~1|obs),# Consider related_temps / related_individuals
                     data=data)
saveRDS(model_taxon, "./Final.Models/Meta_analysis_models/Meta_mods_accl_ratio/meta.acc.taxon.rds")
summary(model_taxon)
# Taxon_plot
orchard_plot(model_taxon, group = "study_ID", mod="class", xlab="ARR", data = data, legend.pos = "bottom.right") +  
  theme(axis.text.y = element_text(angle = 45),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

i2_ml(model_taxon)

# dropping tuatara for figure purposes
Taxon_fig_data <- data %>% filter(class != "tuatara")
fig_model_taxon <- rma.mv(ARR~class-1,  V=Var_ARR,method="REML",test="t",dfs="contain",
                     random=list(~1|study_ID, 
                                 ~1|genus_species,
                                 ~1|obs),# Consider related_temps / related_individuals
                     data=Taxon_fig_data)
Taxon_plot <- orchard_plot(fig_model_taxon, group = "study_ID", mod="class", xlab="ARR", 
                           data = Taxon_fig_data, legend.pos = "bottom.right") +  
  scale_fill_manual(values = c("cornsilk4", "darkgoldenrod1","darkgreen"))+ 
  scale_color_manual(values = c("cornsilk4", "darkgoldenrod1","darkgreen"))+
  theme(axis.text.y = element_text(angle = 45),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

##########
# Figure
##########
# Arranging data for CTMax/Tpref and Overall for final figure
p1 <- mod_results(model_trait, group = "study_ID", mod="1",  data = data)
p2 <- mod_results(model_trait, group = "study_ID", mod="trait", data = data)
p1p2 <- submerge(p1, p2)

# Overall figure
Overall <- orchard_plot(p1p2, group = "study_ID", xlab = "ARR", data = data, 
                        legend.pos = "none") +
  scale_x_discrete(labels = c('Tpref','CTmax','Overall'))+
  scale_fill_manual(values = c("deepskyblue2", "firebrick","grey20"))+ 
  scale_color_manual(values = c("deepskyblue2", "firebrick","grey20"))+
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle = 45))

# Final figure  
# left col
left_col <- plot_grid(Overall,Lifestage_plot,LatLong_plot, labels = c('A', 'B', 'C'), label_size = 12, ncol = 1, align = 'hv', axis = 'l')
right_col <- plot_grid(Taxon_plot, labels = c('D'), label_size = 12)

Fig1 <- plot_grid(left_col, right_col)  
Fig1  

# Figure S2 - metanalysis by species
Fig.S2 <- orchard_plot(model_species, group = "study_ID", mod="genus_species", xlab="AAR", data = data) + 
  theme(axis.text.y = element_text(angle = 35),
        text = element_text(size = 12))  

# Figure S3
model.bias <- rma.mv(yi ~ 0 + scale(temp_diff, scale=FALSE), random = list(~1|study_ID, ~1| species, ~1|obs), V = vi, data = data)
# Funnel plot
funnel( yaxis="seinv", model.bias, legend=FALSE, digit = 2, las = 1))