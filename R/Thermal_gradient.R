########
# Tpref: L. delicata 2020
########
## Packages
pacman::p_load(brms,bayesplot,lme4,tidyverse, latex2exp, performance, plotrix, ggforce, cowplot, emmeans, ggpubr, data.table, rstatix)

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
legend_title <- "Resource Treatment"
# Tpref - bring in data
Tpref.data.raw <- read.csv("./Final.Analysis.Data/Tpref_datasheet_2020.csv") %>%
  select(c("lane_number","bd_liz_id","toeClip","sex",
           "mating_encl_ID",	"body_size", "trt","mean_temp")) %>% 
  mutate(temp = gsub("[AC]-", "", trt),
         treat = gsub("-.*", "", trt)) 
Tpref.data.raw$treat<- as.factor(Tpref.data.raw$treat)
Tpref.data.raw$temp<- as.factor(Tpref.data.raw$temp)
# Tpref- model for plot
m_tpref <- lm(mean_temp ~ scale(body_size) + sex + temp + treat + temp:treat, data = Tpref.data.raw)
m_tpref_emm <- emmeans(m_tpref, specs = c("temp", "treat")) %>% 
  summary() %>%
  data.table()


# Tpref - within treatment comparison for ggplot
Tpref.within.stat.test <- Tpref.data.raw %>% group_by(temp) %>% 
  emmeans_test(mean_temp ~treat, detailed = TRUE) %>% 
  add_xy_position(x = "temp", dodge = 0.8) 
Tpref.within.stat.test<- Tpref.within.stat.test %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  select(temp, estimate, conf.low, conf.high, p)

# Tpref - across treatment comparison for ggplot
Tpref.across.stat.test <- Tpref.data.raw %>% 
  emmeans_test(mean_temp ~temp, detailed = TRUE) %>% 
  add_xy_position(x = "temp", dodge = 0.8) 
Tpref.across.stat.test <- Tpref.across.stat.test %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  select(term, estimate, conf.low, conf.high, p)

# Tpref - figure
pd = position_dodge(.8)
background <- data.frame(lower = seq( 0  , 3), 
                         upper = seq( 1.5, 4.5),
                         col = c("cold", "hot"))
Tpref_fig <-ggplot(m_tpref_emm, aes(x = temp,y= emmean, group = treat, color = treat)) +
  geom_rect(data = m_tpref_emm, aes(x = temp, y=emmean), xmin = as.numeric(m_tpref_emm$temp[[1]]) - 0.6,
            xmax = as.numeric(m_tpref_emm$temp[[1]]) + 0.5,
            ymin = 0, ymax = 46, alpha=0.05, fill="lightblue")+
  geom_rect(data = m_tpref_emm, aes(x = temp, y=emmean), xmin = as.numeric(m_tpref_emm$temp[[2]]) - 0.5,
            xmax = as.numeric(m_tpref_emm$temp[[2]]) + 0.6,
            ymin = 0, ymax = 46, alpha=0.05, fill="brown2")+
  geom_point(shape = 19, size  = 8, position = pd) +
  geom_errorbar(aes(ymin  = lower.CL,
                    ymax  = upper.CL),
                width = 0.2,
                size  = 0.7,
                position = pd)+
  scale_color_manual(legend_title, values=c("black", "goldenrod"), labels=c('Yolk removal', 'Control'))+
  theme_classic()+
  theme(legend.position = "none",
        text=element_text(size=14), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab(bquote(T[Pref]~(degree*C)))+
  scale_y_continuous(breaks=seq(26, 35, 2), limits = c(26,35), expand = expansion(mult = c(0, 0.1))) 


# CTMax (need to bring in data)
CT.data.raw <- read.csv("./Final.Analysis.Data/CTmax_datasheet_2020.csv") %>%
  mutate(temp = gsub("[AC]-", "", trt),
         treat = gsub("-.*", "", trt), 
         age_class = "Juv",
         year = "2020",
         spp = "delicata")

CT.data.raw$treat<- as.factor(CT.data.raw$treat)
CT.data.raw$temp<- as.factor(CT.data.raw$temp)
# CTMax- model for plot
m_CTmax <- lm(CTmax ~ scale(mass) + sex + temp + treat + temp:treat, data = CT.data.raw)
m_CTmax_emm <- emmeans(m_CTmax, specs = c("temp", "treat")) %>% 
  summary() %>%
  data.table()

# CTMax - within treatment comparison for ggplot
CTmax.within.stat.test <- CT.data.raw %>% group_by(temp) %>% 
  emmeans_test(CTmax ~treat, detailed = TRUE) %>% 
  add_xy_position(x = "temp", dodge = 0.8) 
CTmax.within.stat.test <- CTmax.within.stat.test %>% 
  mutate(across(where(is.numeric), round, 2))%>% 
  select(temp, estimate, conf.low, conf.high, p)

# CTMax - across treatment comparison for ggplot
CTmax.across.stat.test <- CT.data.raw %>% 
  emmeans_test(CTmax ~temp, detailed = TRUE) %>% 
  add_xy_position(x = "temp", dodge = 0.8)
CTmax.across.stat.test <- CTmax.across.stat.test%>% 
  mutate(across(where(is.numeric), round, 2))%>% 
  select(term, estimate, conf.low, conf.high, p)

# CTMax Figure
CTmax_fig <- ggplot(m_CTmax_emm, aes(x = temp,y= emmean,color = treat)) +
  geom_rect(data = m_tpref_emm, aes(x = temp, y=emmean), xmin = as.numeric(m_tpref_emm$temp[[1]]) - 0.6,
            xmax = as.numeric(m_tpref_emm$temp[[1]]) + 0.5,
            ymin = 0, ymax = 46, alpha=0.05, fill="lightblue")+
  geom_rect(data = m_tpref_emm, aes(x = temp, y=emmean), xmin = as.numeric(m_tpref_emm$temp[[2]]) - 0.5,
            xmax = as.numeric(m_tpref_emm$temp[[2]]) + 0.6,
            ymin = 0, ymax = 46, alpha=0.05, fill="brown2")+
  geom_point(shape = 19, size  = 8, position = pd) +
  geom_errorbar(aes(ymin  = lower.CL,
                    ymax  = upper.CL),
                width = 0.2,
                size  = 0.7,
                position = pd)+
  scale_color_manual(legend_title, values=c("black", "goldenrod"), labels=c('Yolk removal', 'Control'))+
  theme_classic()+
  theme(legend.position = c(.98, .2),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(5, 5, 5, 5),
        text=element_text(size=14),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), 
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab(bquote(CT[Max]~(degree*C)))+
  scale_y_continuous(breaks=seq(40, 45, 2), limits = c(40,45)) 


# Combining both plots
Temp.Figure <- plot_grid(Tpref_fig, CTmax_fig, labels = c('A', 'B'))