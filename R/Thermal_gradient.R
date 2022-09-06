########
# Tpref: L. delicata 2020
########
## Packages
pacman::p_load(brms,bayesplot,lme4,tidyverse)
## Data
# bring in data and arrange for analysis
Tpref.data.raw <- read.csv("./Final.Analysis.Data/Tpref_datasheet_2020.csv") %>%
  select(c("lane_number","bd_liz_id","toeClip","sex",
           "mating_encl_ID",	"body_size", "trt","mean_temp")) %>% 
  mutate(temp = gsub("[AC]-", "", trt),
         treat = gsub("-.*", "", trt)) 
  
## residual function - used to check residuals for all models
residuals_brms <- function(model, data){
  fitted <- predict(model)[1:dim(data)[1], 1] # Predict mean for each data
  e <- with(data, mean_temp) - fitted
  return(e)
}


####################
# Model 1 & checks(lags, residuals, r2, summary)
####################
model1 <- brms::brm(mean_temp ~ 1 , data = Tpref.data.raw, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000, cores = 4)
model1 <- add_criterion(model1, c("waic", "loo"))
# Checks and plots
m1.draws <- as.array(model1)
mcmc_acf(m1.draws,  pars = c("b_Intercept"), lags =10)
plot(model1)
# Residuals
resid.1 <- residuals_brms(model1, Tpref.data.raw)
hist(resid.1)
# R2 and summary of full model
bayes_R2(model1)
summary(model1)
# save model
saveRDS(model1, "./Final.Models/Tpref_models/Tpref.m1.rds")

####################
# Model 2 & checks(lags, residuals, r2, summary)
####################
model2 <- brms::brm(mean_temp ~ scale(body_size), data = Tpref.data.raw, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000, cores =4)
model2 <- add_criterion(model2, c("waic", "loo"))
# Checks and plots
m2.draws <- as.array(model2)
mcmc_acf(m2.draws,  pars = c("b_Intercept", "b_scalebody_size"), lags =10)
plot(model2)
# Residuals
resid.2 <- residuals_brms(model2, Tpref.data.raw)
hist(resid.2)
# R2 and summary of full model
bayes_R2(model2)
summary(model2)
# save model
saveRDS(model2, "./Final.Models/Tpref_models/Tpref.m2.rds")

####################
# Model 3 & checks(lags, residuals, r2, summary)
####################
model3<- brms::brm(mean_temp ~ scale(body_size) +  temp , data = Tpref.data.raw, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000, core =4)
model3 <- add_criterion(model3, c("waic", "loo"))
# Checks and plots
m3.draws <- as.array(model3)
mcmc_acf(m3.draws,  pars = c("b_Intercept", "b_scalebody_size", "b_temp28"), lags =10)
plot(model3)
# Residuals
resid.3 <- residuals_brms(model3, Tpref.data.raw)
hist(resid.3)
# R2 and summary of full model
bayes_R2(model3)
summary(model3)
# save model
saveRDS(model3, "./Final.Models/Tpref_models/Tpref.m3.rds")

####################
# Model 4 & checks(lags, residuals, r2, summary)
####################
model4 <- brms::brm(mean_temp ~ scale(body_size) +  treat , data = Tpref.data.raw, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000, cores = 4)
model4 <- add_criterion(model4, c("waic", "loo"))
# Checks and plots
m4.draws <- as.array(model4)
mcmc_acf(m4.draws,  pars = c("b_Intercept", "b_scalebody_size", "b_treatC"), lags =10)
plot(model4)
# Residuals
resid.4 <- residuals_brms(model4, Tpref.data.raw)
hist(resid.4)
# R2 and summary of full model
bayes_R2(model4)
summary(model4)
# save model
saveRDS(model4, "./Final.Models/Tpref_models/Tpref.m4.rds")

####################
# Model 5 & checks(lags,residuals, r2, summary)
####################
model5 <- brms::brm(mean_temp ~ scale(body_size) +  temp + treat, data = Tpref.data.raw, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000, cores = 4)
model5 <- add_criterion(model5, c("waic", "loo"))
# Checks and plots
m5.draws <- as.array(model5)
mcmc_acf(m5.draws,  pars = c("b_Intercept", "b_scalebody_size", "b_temp28","b_treatC"), lags =10)
plot(model5)
# Residuals
resid.5 <- residuals_brms(model5, Tpref.data.raw)
hist(resid.5)
# R2 and summary of full model
bayes_R2(model5)
summary(model5)
# save model
saveRDS(model5, "./Final.Models/Tpref_models/Tpref.m5.rds")

####################
# Model 6 & checks(lags, r2, summary)
####################
model6 <- brms::brm(mean_temp ~ scale(body_size) + temp + treat + temp:treat, data = Tpref.data.raw, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000, cores = 4)
model6 <- add_criterion(model6, c("waic", "loo"))
# Checks and plots
m6.draws <- as.array(model6)
mcmc_acf(m6.draws,  pars = c("b_Intercept", "b_scalebody_size", "b_temp28","b_treatC", "b_temp28:treatC"), lags =10)
plot(model6)
# Residuals
resid.6 <- residuals_brms(model5, Tpref.data.raw)
hist(resid.6)
# R2 and summary of full model
bayes_R2(model6)
summary(model6)
# save model
saveRDS(model6, "./Final.Models/Tpref_models/Tpref.m6.rds")

####################
# Model 7 & checks(lags, r2, summary)
####################
model7 <- brms::brm(mean_temp ~ scale(body_size) + sex, data = Tpref.data.raw, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000, cores = 4)
model7 <- add_criterion(model7, c("waic", "loo"))
# Checks and plots
m7.draws <- as.array(model7)
mcmc_acf(m7.draws,  pars = c("b_Intercept", "b_scalebody_size", "b_sexm"), lags =10)
plot(model7)
# Residuals
resid.7 <- residuals_brms(model7, Tpref.data.raw)
hist(resid.7)
# R2 and summary of full model
bayes_R2(model7)
summary(model7)
# save model
saveRDS(model7, "./Final.Models/Tpref_models/Tpref.m7.rds")


####################
# Model comparison
####################
# Both WAIC and LOO model 1 is most parsimonious model followed by model 2 
waic(model1, model2, model3, model4, model5, model6, model7)
model_weights(model1, model2, model3, model4, model5, model6, model7,weights = "waic")
model_weights(model1, model2, model3, model4, model5, model6, model7,weights = "loo")







