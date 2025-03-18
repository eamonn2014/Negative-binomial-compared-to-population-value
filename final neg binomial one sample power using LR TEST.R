# Clear environment and set seed for reproducibility
rm(list=ls())
set.seed(1234)

# Load required packages
library(MASS)    # For glm.nb
library(lmtest)  # For lrtest

# Parameters
true_mean <- 0.58    # Alternative mean (H1)
null_mean <- 1       # Null mean (H0)
r1 <- 1.78           # NB dispersion parameter (size)
alpha <- 0.01        # One-sided significance level
nsim <- 1000         # Number of simulations for stability
fup <- 1             # Follow-up time

# Test Functions
#------------------------------------------------------------------------------------
poisson_test <- function(x, mu0 = 1) {
  test <- poisson.test(sum(x), length(x), r = mu0, alternative = "less")
  return(test$p.value < alpha)
}

#------------------------------------------------------------------------------------
nb_test <- function(x, mu0 = 1, fup = 1) {
  logtime <- rep(log(fup), length(x))
  mod <- tryCatch(
    glm.nb(x ~ 1 + offset(logtime)),
    error = function(e) return(NULL),
    warning = function(w) return(NULL)
  )
  if (is.null(mod)) return(NA)  # Return NA instead of FALSE
  coef_value <- coef(mod)[1]
  se_value <- sqrt(vcov(mod)[1, 1])
  log_mu0 <- log(mu0)
  z_stat <- (coef_value - log_mu0) / se_value
  p <- pnorm(z_stat, lower.tail = TRUE)
  return(p < alpha)
}

#------------------------------------------------------------------------------------
lr_test <- function(x, mu0 = 1, fup = 1) {
  logtime <- rep(log(fup), length(x))
  mod_alt <- tryCatch(
    glm.nb(x ~ 1 + offset(logtime)),
    error = function(e) return(NULL),
    warning = function(w) return(NULL)
  )
  if (is.null(mod_alt)) return(NA)
  theta <- mod_alt$theta  
  ll_null <- sum(dnbinom(x, size = theta, mu = mu0 * fup, log = TRUE))
  ll_alt <- as.numeric(logLik(mod_alt))  # Extract numeric log-likelihood
  lr_stat <- 2 * (ll_alt - ll_null)
  p_value <- pchisq(lr_stat, df = 1, lower.tail = FALSE) / 2
  return(p_value < alpha)
}

#------------------------------------------------------------------------------------
lr_test_fixed <- function(x, mu0 = 1, fup = 1) {
  logtime <- rep(log(fup), length(x))
  mod_alt <- tryCatch(
    glm.nb(x ~ 1 + offset(logtime)),
    error = function(e) return(NULL),
    warning = function(w) return(NULL)
  )
  if (is.null(mod_alt)) return(NA)
  theta_fixed <- 1.78  # Fixed theta (true dispersion)
  ll_null <- sum(dnbinom(x, size = theta_fixed, mu = mu0 * fup, log = TRUE))
  ll_alt <- as.numeric(logLik(mod_alt))  
  lr_stat <- 2 * (ll_alt - ll_null)
  p_value <- pchisq(lr_stat, df = 1, lower.tail = FALSE) / 2
  return(p_value < alpha)
}

#------------------------------------------------------------------------------------
calc_power <- function(n, dist = "poisson") {
  rejections <- replicate(nsim, {
    x <- if (dist == "poisson") {
      rpois(n, lambda = true_mean)
    } else {
      rnbinom(n, size = r1, mu = true_mean)
    }
    
    test_func <- switch(dist,
                        "poisson" = poisson_test,
                        "nb_wald" = nb_test,
                        "nb_lr" = lr_test,
                        "nb_lr_fixed" = lr_test_fixed)
    
    test_func(x)
  })
  
  valid_rejections <- rejections[!is.na(rejections)]  # Remove NA values
  power <- mean(valid_rejections)  # Compute power from valid tests
  
  cat(sprintf("n = %d, %s power = %.3f, NA count = %d\n", 
              n, dist, power, sum(is.na(rejections))))
  return(power)
}

# Test range of n values in jumps of 5
n_range <- seq(40, 100, by = 5)
poisson_powers <- sapply(n_range, calc_power, dist = "poisson")
nb_wald_powers <- sapply(n_range, calc_power, dist = "nb_wald")
nb_lr_powers <- sapply(n_range, calc_power, dist = "nb_lr")
nb_lr_fixed_powers <- sapply(n_range, calc_power, dist = "nb_lr_fixed")

# Find n where power >= 0.8
poisson_n <- ifelse(any(poisson_powers >= 0.8), n_range[min(which(poisson_powers >= 0.8))], NA)
nb_wald_n <- ifelse(any(nb_wald_powers >= 0.8), n_range[min(which(nb_wald_powers >= 0.8))], NA)
nb_lr_n <- ifelse(any(nb_lr_powers >= 0.8), n_range[min(which(nb_lr_powers >= 0.8))], NA)
nb_lr_fixed_n <- ifelse(any(nb_lr_fixed_powers >= 0.8), n_range[min(which(nb_lr_fixed_powers >= 0.8))], NA)

# Results
cat("Poisson: n =", poisson_n, "Power =", round(poisson_powers[which(n_range == poisson_n)], 3), "\n")
cat("Neg Bin Wald (r = 1.78): n =", nb_wald_n, "Power =", round(nb_wald_powers[which(n_range == nb_wald_n)], 3), "\n")
cat("Neg Bin LR (r = 1.78): n =", nb_lr_n, "Power =", round(nb_lr_powers[which(n_range == nb_lr_n)], 3), "\n")
cat("Neg Bin LR Fixed θ (r = 1.78): n =", nb_lr_fixed_n, "Power =", 
    round(nb_lr_fixed_powers[which(n_range == nb_lr_fixed_n)], 3), "\n")

# Plot power curves with fourth curve (purple)
plot(n_range, poisson_powers, type = "l", col = "green", ylim = c(0, 1),
     xlab = "Sample Size (n)", ylab = "Empirical Power",
     main = "Power vs. Sample Size (alpha = 0.01, μ = 0.58)")
lines(n_range, nb_wald_powers, col = "blue")
lines(n_range, nb_lr_powers, col = "red")
lines(n_range, nb_lr_fixed_powers, col = "purple")
abline(h = 0.8, lty = 2)
legend("bottomright", 
       legend = c("Poisson (λ = 0.58)", "NB Wald (r = 1.78)", 
                  "NB LR (r = 1.78)", "NB LR Fixed θ (r = 1.78)"),
       col = c("green", "blue", "red", "purple"), lty = 1)
