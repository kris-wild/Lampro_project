########
# Meta-analysis
########
## Packages
# Load libraries
pacman::p_load(devtools, dplyr, metafor, patchwork, R.rsp, emmeans, ggdist, rotl, ape, orchaRd, patchwork, latex2exp, cowplot, ggplot2)

## data
# bring in data and arrange for analysis
data <- read.csv("./Final.Analysis.Data/Final_Meta_results_extracted_data.csv")

### ### ### 
### Meta-analysis with acclimation response ratio #### 
### ### ### 

data_cohort <- group_by(data,cohort_ID) %>% summarise(n_cohort=n()) # Identify repeated observations on the same individuals
data<- left_join(data, data_cohort) # Join with main dataset


data <- data %>%  mutate(ARR = ((trt_mean-c_mean)/(trt_temp-c_temp)), # Calculate ARR
                   
                   Var_ARR = ifelse(n_cohort=="1", # If independent, use this formula
                                     ((1/(trt_temp - c_temp))^2*(c_sd^2/c_n+trt_sd^2/trt_n)), # use this formula 
                                       # If non-independent, use this formula
                                            ((1/(trt_temp - c_temp))^2*((c_sd^2 + trt_sd^2 - 2*0.5*c_sd*trt_sd)/(c_n + trt_n)))), 
                   
                   precision =(1/sqrt(Var_ARR))) # Calculate precision (1/SE)


# acclimation response ratio figure
ggplot(data) + 
  geom_vline(xintercept=0, alpha=0.5, linetype=2)+
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

## Create VCV matrix for sampling errors because the same data was used to calculate multiple effect sizes (stepwise comparisons; shared_trt_ID)
names(data)<-make.names(names(data))
data<-as.data.frame(data)
VCV_ARR <- metaAidR::make_VCV_matrix(data, V = "Var_ARR", cluster = "shared_trt_ID", obs = "obs_ID")

## save final "data" df for rmd results  
data_final <- saveRDS(data, "./Final.Analysis.Data/Meta_df_rmdfile.rds")

###################################################  
##### Phylo, spp, observation - intercept mods
###################################################  
# 1) Model ALL with phylo: Both Tpref and CTmax analysed together, with all random effects
data$spp <- data$genus_species

model_all_rand<-rma.mv(ARR~1, 
                       V=VCV_ARR,# Consider using the VCV matrix of ARR with correlated errors
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|genus_species, #species name with covariance matrix
                                   ~1|spp, # species name
                                   ~1|obs_ID),
                       R=list(genus_species=phylo_matrix),
                       data=data)

# check phylogeny  
orchaRd::i2_ml(model_all_rand) #  phylogeny explains virtually no variance, so we can remove it
I2_all_phylo <- as.data.frame(round(i2_ml(model_all_rand),digits = 2))
I2_all_phylo
I2_all_phylo_predict <- as.data.frame(predict(model_all_rand))
I2_all_phylo_predict
AIC.rma(model_all_rand)
saveRDS(model_all_rand, "./Final.Models/Meta_analysis_models/meta.acc.phylo.allrand.rds")

# 2) model removing species name as random variable but keeping phylo matrix
model_all_no_spp<-rma.mv(ARR~1, 
                       V=VCV_ARR,
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|genus_species,
                                   ~1|obs_ID),
                       R=list(genus_species=phylo_matrix),
                       data=data)
# check if spp drives varience
orchaRd::i2_ml(model_all_no_spp) #  phylogeny explains virtually no variance, so we can remove it
I2_all_phylo <- as.data.frame(round(i2_ml(model_all_no_spp),digits = 2))
I2_all_phylo_predict <- as.data.frame(predict(model_all_no_spp))
I2_all_phylo_predict
AIC.rma(model_all_no_spp)
saveRDS(model_all_no_spp, "./Final.Models/Meta_analysis_models/meta.acc.phylo.nospp.rds")

# 3) Drop phylogeny and just check if species is important. Basically the same so keep species in.
model_all_rand_spp<-rma.mv(ARR~1, 
                       V=VCV_ARR,
                       method="REML",
                       test="t",
                       dfs="contain",
                       random=list(~1|study_ID, 
                                   ~1|genus_species,
                                   ~1|obs_ID),
                       data=data)
summary(model_all_rand_spp)     
i2_ml(model_all_rand_spp) 
AIC.rma(model_all_rand_spp)
I2_model_all_rand_spp_predict <- as.data.frame(predict(model_all_rand_spp))

#############################################
# Fig2A) Intercept model: **without phylogeny
#############################################
int_model<-rma.mv(ARR~1, 
                  V=VCV_ARR,
                  method="REML",
                  test="t",
                  dfs="contain",
                  random=list(~1|study_ID, 
                              ~1|genus_species,
                              ~1|obs_ID),
                  data=data)


saveRDS(int_model, "./Final.Models/Meta_analysis_models/meta.acc.no.phylo.rds")
# model checks Intercept model
summary(int_model)     
i2_ml(int_model) # Lots of heterogeneity
predict(int_model) # Prediction intervals
AIC.rma(int_model)

###############################
# Fig2A): Trait differences mod
###############################
model_trait<- rma.mv(ARR~trait-1, 
                     V=VCV_ARR,
                     method="REML",
                     test="t",
                     dfs="contain",
                     random=list(~1|study_ID,
                                 ~1|genus_species,
                                 ~1|obs_ID),
                     data=data)

## Lets use RVE to account for non-independence
model_trait_RVE <- robust(model_trait, cluster = data$cluster)

saveRDS(model_trait, "./Final.Models/Meta_analysis_models/meta.acc.trait.rds")
summary(model_trait)
orchard_plot(model_trait, group = "study_ID", mod="trait", xlab="ARR", data = data)
i2_ml(model_trait)
I2_model_trait <- as.data.frame(round(i2_ml(model_trait),digits = 2))
I2_model_trait
mod_results(model_trait, mod="trait", group="study_ID", data=data) # Prediction intervals

#########################
## Fig2B) Life stage mod
#########################
data$age <- factor(data$age, levels=c("adult", "juvenile", "hatchling"))
model_age<- rma.mv(ARR~age-1, 
                   V=VCV_ARR,
                   method="REML",
                   test="t",
                   dfs="contain",
                   random=list(~1|study_ID, 
                               ~1|genus_species,
                               ~1|obs_ID),
                   data=data)

saveRDS(model_age, "./Final.Models/Meta_analysis_models/meta.acc.age.rds")
summary(model_age)
Lifestage_plot <- orchard_plot(model_age, group = "study_ID", mod="age", xlab="ARR", data = data, legend.pos = "none", trunk.size=6) + 
  scale_fill_manual(values = c("deepskyblue2","firebrick", "darkgoldenrod1"))+ 
  scale_color_manual(values = c("deepskyblue2","firebrick", "darkgoldenrod1"))+
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle = 45))
i2_ml(model_age)
mod_results(model_age, mod="age", group="study_ID", data=data) # Prediction intervals

####################
# Fig2C) Temperate vs Tropical mod
####################
model_zone<- rma.mv(ARR~zone-1, 
                    V=VCV_ARR,
                    method="REML",
                    test="t",
                    dfs="contain",
                    random=list(~1|study_ID, 
                                ~1|genus_species,
                                ~1|obs_ID),
                    data=data)

saveRDS(model_zone, "./Final.Models/Meta_analysis_models/meta.acc.zone.rds")
summary(model_zone)
LatLong_plot <-orchard_plot(model_zone, group = "study_ID", mod="zone", xlab="ARR", data = data, legend.pos = "none", trunk.size=6) +
  scale_fill_manual(values = c("deepskyblue2","firebrick"))+ 
  scale_color_manual(values = c("deepskyblue2","firebrick"))+
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text.y = element_text(angle = 45))
i2_ml(model_zone)
mod_results(model_zone, mod="zone", group="study_ID", data=data) # Prediction intervals

####################
## Fig2D) major taxonomic group mod
####################
model_taxon<- rma.mv(ARR~class-1, 
                     V=VCV_ARR,# Consider using the VCV matrix of ARR with correlated errors
                     method="REML",
                     test="t",
                     dfs="contain",
                     random=list(~1|study_ID, 
                                 ~1|genus_species,
                                 ~1|obs_ID),
                     data=data)

saveRDS(model_taxon, "./Final.Models/Meta_analysis_models/meta.acc.taxon.rds")
summary(model_taxon)
mod_results(model_taxon, mod="class", group="study_ID", data=data) # Prediction intervals

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
VCV_ARRt <- metaAidR::make_VCV_matrix(Taxon_fig_data, V = "Var_ARR", cluster = "shared_trt_ID", obs = "obs_ID")
fig_model_taxon <- rma.mv(ARR~class-1,  
                          V=VCV_ARRt,
                          method="REML",
                          test="t",
                          dfs="contain",
                     random=list(~1|study_ID, 
                                 ~1|genus_species,
                                 ~1|obs_ID),
                     data=Taxon_fig_data)

Taxon_plot <- orchard_plot(fig_model_taxon, group = "study_ID", mod="class", xlab="ARR", 
                           data = Taxon_fig_data, legend.pos = "bottom.right", trunk.size=6) +  
  scale_fill_manual(values = c("cornsilk4", "darkgoldenrod1","darkgreen"))+ 
  scale_color_manual(values = c("cornsilk4", "darkgoldenrod1","darkgreen"))+
  theme(axis.text.y = element_text(angle = 45),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


####################
# Mod bias Check
####################

model.bias1   <- rma.mv(ARR~ sqrt(Var_ARR),  # Standard error
                          V=VCV_ARR,
                          method="REML",
                          test="t",
                          dfs="contain",
                          random=list(~1|study_ID, 
                                      ~1|genus_species,
                                      ~1|obs_ID),
                          data=data)

summary(model.bias1) 

model.bias2   <- rma.mv(ARR~ Var_ARR,  # Variance
                        V=VCV_ARR,
                        method="REML",
                        test="t",
                        dfs="contain",
                        random=list(~1|study_ID, 
                                    ~1|genus_species,
                                    ~1|obs_ID),
                        data=data)
summary(model.bias2)

# Funnel plot
funnel( yaxis="seinv", model.bias1, legend=FALSE, digit = 2, las = 1)
saveRDS(model.bias1, "./Final.Models/Meta_analysis_models/mod_bias.rds")






##########
# Figure2 - FINAL
##########
# Arranging data for CTMax/Tpref and Overall for final figure
p1 <- mod_results(model_trait, group = "study_ID", mod="1",  data = data)
p2 <- mod_results(model_trait, group = "study_ID", mod="trait", data = data)
p1p2 <- submerge(p1, p2)

# Overall figure
Overall <- orchard_plot(p1p2, group = "study_ID", xlab = "ARR", data = data, 
                        legend.pos = "none", trunk.size=6) +
  scale_x_discrete(labels = c('Tpref','CTmax','Overall'))+
  scale_fill_manual(values = c("deepskyblue2", "firebrick","grey20"))+ 
  scale_color_manual(values = c("deepskyblue2", "firebrick","grey20"))+
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.y = element_text(angle = 45))

# Final figure  
# left col
left_col <- plot_grid(Overall,
                      Lifestage_plot,
                      LatLong_plot, 
                      labels = c('A', 'B', 'C'), 
                      label_size = 12, 
                      ncol = 1, 
                      align = 'hv', 
                      axis = 'l')
right_col <- plot_grid(Taxon_plot, labels = c('D'), label_size = 12)

Fig2 <- plot_grid(left_col, right_col)  
Fig2  

