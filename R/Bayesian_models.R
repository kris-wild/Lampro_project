########
# Bayesian Models
########
# Packages
pacman::p_load(brms)

BM <- seq(5, 30, length.out = 50)
int <- 5
B_mass <- 0.10
B_trt <- 0.8
trt <- rep(c(0,1), each = 25)
CTMax = int + B_mass*BM + B_trt*trt + rnorm(50, 0, 1)
plot(CTMax ~ BM)

model <- glm(CTMax ~ BM + trt)
summary(model)



data <- data.frame(CTMax = CTMax, BM = BM, trt = trt)
model <- brms::brm(CTMax ~ BM + trt, data = data, family = gaussian(), chains = 4, iter = 2000, thin = 1, warmup = 1000)
summary(model)

plot(model)
posterior <- posterior_samples(model)
head(posterior)

b_0 <- posterior[,"b_Intercept"]
head(b_0)
b_1 <- posterior[,"b_Intercept"] + posterior[,"b_trt"]

contrast <- b_1-b_0
mean(contrast)
p_value <- sum(table(contrast) <= 0) / (length(contrast)-1)
posterior_bm <- posterior_samples(model)[,"b_BM"]
hist(posterior_bm)

mean(posterior_bm)
sd(posterior_bm)
quantile(posterior_bm, c(0.30,0.80))

