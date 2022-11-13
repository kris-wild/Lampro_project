########
# CT Max L. delicata juviniles 2020
########
## Packages
pacman::p_load(brms,bayesplot,lme4,tidyverse,performance, latex2exp)

## Data - bring in data and arrange for analysis
# 2020 - juvenile L. delicata
CT.data.raw <- read.csv("./Final.Analysis.Data/CTmax_datasheet_2020.csv") %>%
  mutate(temp = gsub("[AC]-", "", trt),
         treat = gsub("-.*", "", trt), 
         age_class = "Juv",
         year = "2020",
         spp = "delicata")


####################
# Model 1 
####################
model1 <- lm(CTmax ~ 1, data = CT.data.raw)
m1.summary <- summary(model1) 
# save model
saveRDS(model1, "./Final.Models/CTmax_models/CT.Max.m1.rds")


####################
# Model 2 & checks
####################
model2 <- lm(CTmax ~ sex, data = CT.data.raw)
# Checks and plots
m2.summary <- summary(model2) 
check_model(model2)
# save model
saveRDS(model2, "./Final.Models/CTmax_models/CT.Max.m2.rds")

        
####################
# Model 3 & checks
####################
model3 <- lm(CTmax ~ scale(mass) + sex, data = CT.data.raw)
# Checks and plots
m3.summary <- summary(model3) 
check_model(model3)
# save model
saveRDS(model3, "./Final.Models/CTmax_models/CT.Max.m3.rds")
 

####################
# Model 4 & checks
####################
model4 <- lm(CTmax ~ scale(mass) + sex + temp, data = CT.data.raw)
# Checks and plots
m4.summary <- summary(model4) 
check_model(model4)
# save model
saveRDS(model4, "./Final.Models/CTmax_models/CT.Max.m4.rds")


####################
# Model 5 & checks
####################
model5 <- lm(CTmax ~ scale(mass) + sex + treat, data = CT.data.raw)
# Checks and plots
m5.summary <- summary(model5) 
check_model(model5)
# save model
saveRDS(model5, "./Final.Models/CTmax_models/CT.Max.m5.rds")


####################
# Model 6 & checks
####################
model6 <- lm(CTmax ~ scale(mass) + sex + temp + treat + temp:treat, data = CT.data.raw)
# Checks and plots
m6.summary <- summary(model6) 
check_model(model6)
# save model
saveRDS(model6, "./Final.Models/CTmax_models/CT.Max.m6.rds")


