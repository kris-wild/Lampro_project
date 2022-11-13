########
# Tpref: L. delicata 2020
########
## Packages
pacman::p_load(brms,bayesplot,lme4,tidyverse, latex2exp, performance, ggforce, cowplot)

## Data
# bring in data and arrange for analysis
Tpref.data.raw <- read.csv("./Final.Analysis.Data/Tpref_datasheet_2020.csv") %>%
  select(c("lane_number","bd_liz_id","toeClip","sex",
           "mating_encl_ID",	"body_size", "trt","mean_temp")) %>% 
  mutate(temp = gsub("[AC]-", "", trt),
         treat = gsub("-.*", "", trt)) 
  

####################
# Model 1 
####################
model1 <- lm(mean_temp ~ 1 , data = Tpref.data.raw)
m1.summary <- summary(model1) 
saveRDS(model1, "./Final.Models/Tpref_models/Tpref.m1.rds")


####################
# Model 2 & checks
####################
model2 <- lm(mean_temp ~ sex, data = Tpref.data.raw)
# Checks and plots
m2.summary <- summary(model2) 
check_model(model2)
# save model
saveRDS(model2, "./Final.Models/Tpref_models/Tpref.m2.rds")

####################
# Model 3 & checks
####################
model3 <- lm(mean_temp ~ scale(body_size) + sex, data = Tpref.data.raw)
# Checks and plots
m3.summary <- summary(model3) 
check_model(model3)
# save model
saveRDS(model3, "./Final.Models/Tpref_models/Tpref.m3.rds")

####################
# Model 4 & checks
####################
model4 <- lm(mean_temp ~ scale(body_size) + sex +  temp, data = Tpref.data.raw)
# Checks and plots
m4.summary <- summary(model4) 
check_model(model4)
# save model
saveRDS(model4, "./Final.Models/Tpref_models/Tpref.m4.rds")

####################
# Model 5 & checks
####################
model5 <- lm(mean_temp ~ scale(body_size) + sex +  treat, data = Tpref.data.raw)
# Checks and plots
m5.summary <- summary(model5) 
check_model(model5)
# save model
saveRDS(model5, "./Final.Models/Tpref_models/Tpref.m5.rds")

####################
# Model 6 & checks
####################
model6 <- lm(mean_temp ~ scale(body_size) + sex + temp + treat + temp:treat, data = Tpref.data.raw)
# Checks and plots
m6.summary <- summary(model6) 
check_model(model6)
# save model
saveRDS(model6, "./Final.Models/Tpref_models/Tpref.m6.rds")


########
# Thermal Figures (Tpref & CTmax)
########
# Tpref fig
legend_title <- "Resource Treatment"
tpref.fig <- ggplot(Tpref.data.raw, aes(x = temp, y = mean_temp, fill = treat)) +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_violin(trim = FALSE) +
  geom_sina(alpha=0.65)+
  scale_fill_manual(values=c("brown2", "gray82"))+
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size=14))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab(bquote(T[Pref]~(degree*C)))+
  scale_y_continuous(breaks=seq(18, 38, 2), limits = c(18,38)) 


# CTMax (need to bring in data)
CT.data.raw <- read.csv("./Final.Analysis.Data/CTmax_datasheet_2020.csv") %>%
  mutate(temp = gsub("[AC]-", "", trt),
         treat = gsub("-.*", "", trt), 
         age_class = "Juv",
         year = "2020",
         spp = "delicata")
# figure
ctmax.fig <- ggplot(CT.data.raw, aes(x = temp, y = CTmax, fill = treat)) +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_violin(trim = FALSE) +
  geom_sina(alpha=0.65)+
  scale_fill_manual(legend_title, values=c("brown2", "gray82"), labels=c('Yolk ablation', 'Control'))+ 
  theme_classic() +
  theme(legend.position = c(.99, .2),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        text=element_text(size=14))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab(bquote(CT[Max]~(degree*C)))+
  scale_y_continuous(breaks=seq(36, 48, 2), limits = c(36,48)) 


# Combining both plots
Temp.Figure <- plot_grid(tpref.fig, ctmax.fig, labels = c('A', 'B'))
# draw_label("Treatment", x=0.5, y=  0, vjust=-0.5, angle= 0)

