


 #1
 ## too simplisiticc!
 
 bci <- function(n,r, dispersion, fup) {
 
   # n=80
   # r = 0.58  # 
   # fup=1
   # dispersion = 1.78 # k
   
   y <- rnbinom(n, mu = r * fup, size = dispersion)
   
   # table(y)/n
   # barplot(table(y)/n)
 
   mod <- (glm.nb(y ~ 1 )) #, # offset(1))),
   coef_value <- coef(mod)[1]
   se_value <- sqrt(vcov(mod)[1, 1])
   u <- exp (coef_value + se_value * qnorm(.99))[[1]]
              
 return(u)
 
 }
 
 res <- bci(n=80,r=0.58, dispersion=1.78, fup=1)
 
 mean(replicate(1000, bci(n=70,r=0.58, dispersion=1.78, fup=1))<1)
 
 ##---------------------------------------------------------------------------------------------
 
  #2
 library(MASS) # For negative binomial GLM
 library(lmtest) # For likelihood ratio test
 
 # Function to calculate the upper bound of a confidence interval for a negative binomial mean.
 bci <- function(n, r, dispersion, fup) {
   y <- rnbinom(n, mu = r * fup, size = dispersion)
   mod <- glm.nb(y ~ 1)
   coef_value <- coef(mod)[1]
   se_value <- sqrt(vcov(mod)[1, 1])
   u <- exp(coef_value + se_value * qnorm(0.99))[[1]]
   return(u)
 }
 
 # Function to perform a likelihood ratio test for a one-sample negative binomial.
 lr_test <- function(n, r, dispersion, fup, pop_rate) {
   y <- rnbinom(n, mu = r * fup, size = dispersion)
   mod_alt <- glm.nb(y ~ 1)
   mod_null <- glm.nb(y ~ 1, init.theta = dispersion, mu = rep(pop_rate, length(y))) # Corrected line
   lr_result <- lrtest(mod_alt, mod_null)
   p_value <- lr_result$Pr[2]
   return(p_value)
 }
 
 # Power simulation using confidence interval approach
 power_ci <- function(n, r, dispersion, fup, pop_rate, alpha = 0.01, reps = 1000) {
   ci_test_results <- replicate(reps, bci(n, r, dispersion, fup) < pop_rate)
   power_estimate <- mean(ci_test_results)
   return(power_estimate)
 }
 
 # Power simulation using likelihood ratio test
 power_lr <- function(n, r, dispersion, fup, pop_rate, alpha = 0.01, reps = 1000) {
   lr_test_results <- replicate(reps, lr_test(n, r, dispersion, fup, pop_rate) < alpha)
   power_estimate <- mean(lr_test_results)
   return(power_estimate)
 }
 
 # Example usage and power calculation:
 n <- 90
 r <- 0.58
 dispersion <- 1.78
 fup <- 1
 pop_rate <- 1
 alpha <- 0.01
 
 # Power estimation using confidence interval
 power_estimate_ci <- power_ci(n, r, dispersion, fup, pop_rate, alpha)
 print(paste("Power (CI method):", power_estimate_ci))
 
 # Power estimation using likelihood ratio test
 power_estimate_lr <- power_lr(n, r, dispersion, fup, pop_rate, alpha)
 print(paste("Power (LR test method):", power_estimate_lr))
 
 # Example of a single lrtest
 single_lr_test <- lr_test(n=70, r=0.58, dispersion=1.78, fup=1, pop_rate =1)
 print(paste("Single LR test p value:", single_lr_test))
 
 # Example of a single confidence interval test
 single_ci_test <- bci(n=70, r=0.58, dispersion=1.78, fup=1) < 1
 print(paste("single CI test result:", single_ci_test))
 
 
 #---------------------------------------------------------------------------------------------------------------
 #3
 # bootstrap1
 # Load necessary libraries
library(MASS) # For negative binomial GLM
library(lmtest) # For likelihood ratio test

# Function to perform likelihood ratio test for a one-sample negative binomial.
lr_test <- function(n, r, dispersion, fup, pop_rate) {
  # Generate response variable from Negative Binomial distribution
  y <- rnbinom(n, mu = r * fup, size = dispersion)
  
  # Ensure pop_rate is a vector of the same length as y
  pop_rate_vector <- rep(pop_rate, length(y))
  
  # Fit the unconstrained model
  mod_alt <- glm.nb(y ~ 1)
  
  # Fit the constrained model with offset
  mod_null <- glm.nb(y ~ offset(log(pop_rate_vector)), init.theta = dispersion)
  
  # Perform likelihood ratio test between models
  lr_result <- lrtest(mod_alt, mod_null)
  
  # Extract p-value from the test result
  p_value <- lr_result$Pr[2]
  
  return(p_value)
}

# Power simulation using likelihood ratio test
power_lr <- function(n, r, dispersion, fup, pop_rate, alpha = alpha, reps = reps.) {
  # Simulate likelihood ratio tests and calculate power
  lr_test_results <- replicate(reps, lr_test(n, r, dispersion, fup, pop_rate) < alpha)
  power_estimate <- mean(lr_test_results)
  
  return(power_estimate)
}

# Bootstrap-based power estimation
bootstrap_power_lr <- function(n, r, dispersion, fup, pop_rate, alpha = 0.01, reps = reps., 
                               bootstrap_reps = bootstrap_reps.) {
  # Vector to store power estimates from each bootstrap sample
  bootstrap_power_estimates <- numeric(bootstrap_reps)
  
  for (i in 1:bootstrap_reps) {
    # Generate a bootstrap sample of the data
    y_bootstrap <- rnbinom(n, mu = r * fup, size = dispersion)
    
    # Estimate power using the bootstrap sample
    bootstrap_power_estimates[i] <- power_lr(n, r, dispersion, fup, pop_rate, alpha, reps)
  }
  
  # Calculate the mean and confidence interval of the bootstrap power estimates
  mean_power <- mean(bootstrap_power_estimates)
  ci_power <- quantile(bootstrap_power_estimates, c(0.025, 0.975))
  
  return(list(mean_power = mean_power, ci_power = ci_power))
}

# Example usage and power calculation
n <- 90
r <- 0.58
dispersion <- 1.78
fup <- 1
pop_rate <- 1
alpha <- 0.01
reps. = 100
bootstrap_reps. = 100

# Power estimation using likelihood ratio test
power_estimate_lr <- power_lr(n, r, dispersion, fup, pop_rate, alpha)
print(paste("Power (LR test method):", power_estimate_lr))

# Bootstrap-based power estimation
bootstrap_power_estimate <- bootstrap_power_lr(n, r, dispersion, fup, pop_rate, alpha)
print(paste("Bootstrap Power Estimate (mean):", bootstrap_power_estimate$mean_power))
print(paste("Bootstrap Power Estimate (95% CI):", bootstrap_power_estimate$ci_power))
 
 
 
 #------------------------------------------------------------------
# bootstrap 4

library(MASS)
library(lmtest)

# Function to perform likelihood ratio test for a one-sample negative binomial.
lr_test <- function(n, r, dispersion, fup, pop_rate) {
  y <- rnbinom(n, mu = r * fup, size = dispersion)
  pop_rate_vector <- rep(pop_rate, length(y))
  mod_alt <- glm.nb(y ~ 1)
  mod_null <- glm.nb(y ~ offset(log(pop_rate_vector)), init.theta = dispersion)
  lr_result <- lrtest(mod_alt, mod_null)
  p_value <- lr_result$Pr[2]
  return(p_value)
}

# Power simulation using likelihood ratio test (sequential)
power_lr <- function(params) {
  lr_test_results <- replicate(params$reps, lr_test(params$n, params$r, params$dispersion, params$fup, 
                                                    params$pop_rate) < params$alpha)
  power_estimate <- mean(lr_test_results)
  return(power_estimate)
}

# Bootstrap-based power estimation (sequential)
bootstrap_power_lr <- function(params) {
  bootstrap_power_estimates <- replicate(params$bootstrap_reps, {
    y_bootstrap <- rnbinom(params$n, mu = params$r * params$fup, size = params$dispersion)
    power_lr(params)
  })
  mean_power <- mean(bootstrap_power_estimates)
  ci_power <- quantile(bootstrap_power_estimates, c(0.025, 0.975))
  return(list(mean_power = mean_power, ci_power = ci_power))
}

# Parameter object
params <- list(
  n = 90,
  r = 0.58,
  dispersion = 1.78,
  fup = 1,
  pop_rate = 1,
  alpha = 0.01,
  reps = 100,
  bootstrap_reps = 100
)

# Power estimation using likelihood ratio test
power_estimate_lr <- power_lr(params)
print(paste("Power (LR test method):", power_estimate_lr))

# Bootstrap-based power estimation
bootstrap_power_estimate <- bootstrap_power_lr(params)
print(paste("Bootstrap Power Estimate (mean):", bootstrap_power_estimate$mean_power))
print(paste("Bootstrap Power Estimate (95% CI):", bootstrap_power_estimate$ci_power))

#------------------------------------------------------------------------------------------------------

# bootstrap 5

library(MASS)
library(lmtest)

# Function to perform likelihood ratio test for a one-sample negative binomial.
lr_test <- function(n, r, dispersion, fup, pop_rate) {
  y <- rnbinom(n, mu = r * fup, size = dispersion)
  pop_rate_vector <- rep(pop_rate, length(y))
  mod_alt <- glm.nb(y ~ 1)
  mod_null <- glm.nb(y ~ offset(log(pop_rate_vector)), init.theta = dispersion)
  lr_result <- lrtest(mod_alt, mod_null)
  p_value <- lr_result$Pr[2]
  return(p_value)
}

# Power simulation using likelihood ratio test (sequential)
power_lr <- function(params) {
  lr_test_results <- replicate(params$reps, lr_test(params$n, params$r, params$dispersion, 
                                                    params$fup, params$pop_rate) < params$alpha)
  power_estimate <- mean(lr_test_results)
  return(power_estimate)
}

# Bootstrap-based power estimation (sequential)
bootstrap_power_lr <- function(params) {
  bootstrap_power_estimates <- replicate(params$bootstrap_reps, {
    y_bootstrap <- rnbinom(params$n, mu = params$r * params$fup, size = params$dispersion)
    # Update params with the bootstrap sample
    params_bootstrap <- params
    params_bootstrap$y <- y_bootstrap
    power_lr(params_bootstrap)
  })
  mean_power <- mean(bootstrap_power_estimates)
  ci_power <- quantile(bootstrap_power_estimates, c(0.025, 0.975))
  return(list(mean_power = mean_power, ci_power = ci_power))
}

# Parameter object
params <- list(
  n = 90,                # Sample size
  r = 0.58,              # Rate ratio
  dispersion = 1.78,     # Dispersion parameter
  fup = 1,               # Follow-up time
  pop_rate = 1,          # Population rate
  alpha = 0.01,          # Significance level
  reps = 100,            # Number of replicates for power estimation
  bootstrap_reps = 100   # Number of bootstrap replicates
)

# Power estimation using likelihood ratio test
power_estimate_lr <- power_lr(params)
print(paste("Power (LR test method):", power_estimate_lr))

# Bootstrap-based power estimation
bootstrap_power_estimate <- bootstrap_power_lr(params)
print(paste("Bootstrap Power Estimate (mean):", bootstrap_power_estimate$mean_power))
print(paste("Bootstrap Power Estimate (95% CI):", bootstrap_power_estimate$ci_power))



 #--------------------------------------------------------------------------------------------------------------
#  code 6
 # Load necessary libraries
 rm(list=ls())
 library(MASS) # For negative binomial GLM
 library(lmtest) # For likelihood ratio test
 
 # Function to calculate the upper bound of a confidence interval for a negative binomial mean.
 bci <- function(n, r, dispersion, fup) {
   # Generate response variable from Negative Binomial distribution
   y <- rnbinom(n, mu = r * fup, size = dispersion)
   
   # Fit Negative Binomial model
   mod <- glm.nb(y ~ 1)
   
   # Extract coefficient and standard error
   coef_value <- coef(mod)[1]
   se_value <- sqrt(vcov(mod)[1, 1])
   
   # Calculate upper bound of confidence interval
   u <- exp(coef_value + se_value * qnorm(0.99))
   
   return(u)
 }
 
 # Function to perform likelihood ratio test for a one-sample negative binomial.
 lr_test <- function(n, r, dispersion, fup, pop_rate) {
   # Generate response variable from Negative Binomial distribution
   y <- rnbinom(n, mu = r * fup, size = dispersion)
   
   # Ensure pop_rate is a vector of the same length as y
   pop_rate_vector <- rep(pop_rate, length(y))
   
   # Fit the unconstrained model
   mod_alt <- glm.nb(y ~ 1)
   
   # Fit the constrained model with offset
   mod_null <- glm.nb(y ~ offset(log(pop_rate_vector)), init.theta = dispersion)
   
   # Perform likelihood ratio test between models
   lr_result <- lrtest(mod_alt, mod_null)
   
   # Extract p-value from the test result
   p_value <- lr_result$Pr[2] # 1 or 0?
   
   return(p_value)
 }
 
 # Power simulation using confidence interval approach
 power_ci <- function(n, r, dispersion, fup, pop_rate, alpha = 0.01, reps = 1000) {
   # Simulate confidence intervals and check if upper bound is less than pop_rate
   ci_test_results <- replicate(reps, bci(n, r, dispersion, fup) < pop_rate)
   power_estimate <- mean(ci_test_results)
   
   return(power_estimate)
 }
 
 # Power simulation using likelihood ratio test
 power_lr <- function(n, r, dispersion, fup, pop_rate, alpha = 0.01, reps = 1000) {
   # Simulate likelihood ratio tests and calculate power
   lr_test_results <- replicate(reps, lr_test(n, r, dispersion, fup, pop_rate) < alpha)
   power_estimate <- mean(lr_test_results)
   
   return(power_estimate)
 }
 
 # Example usage and power calculation
 n <- 90
 r <- 0.58
 dispersion <- 1.78
 fup <- 1
 pop_rate <- 1
 alpha <- 0.01
 
 # Power estimation using confidence interval
 power_estimate_ci <- power_ci(n, r, dispersion, fup, pop_rate, alpha)
 print(paste("Power (CI method):", power_estimate_ci))
 
 # Power estimation using likelihood ratio test
 power_estimate_lr <- power_lr(n, r, dispersion, fup, pop_rate, alpha)
 print(paste("Power (LR test method):", power_estimate_lr))
 
 #---------------------------------------------------------------------------------------------
 
 
 #7
 
 ###the LR test is superior to this:::
 
 
 bci <- function(n,r, dispersion, fup) {
   
   y <- rnbinom(n, mu = r * fup, size = dispersion)
   mod <- (glm.nb(y ~ 1 )) #, # offset(1))),
   u <- confint(mod, level = 0.99)[2][[1]]
   return(u)
   
 }
 
 res <- bci(n=80,r=0.58, dispersion=1.78, fup=1)
 
 mean(replicate(1000, bci(n=78, r=0.58, dispersion=1.78, fup=1))<0)
 
 
 
 
 
 
 
 
  