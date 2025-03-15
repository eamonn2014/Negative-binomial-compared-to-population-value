



# Parameters
mean <- 0.58
r1 <- 1.78          # First dispersion parameter
#r2 <- 1 / 1.78      # Second dispersion parameter ≈ 0.562

par(mfrow = c(1, 1))
# Define x range
x <- 0:9

# Calculate PMF for all distributions
pmf_nb1 <- dnbinom(x, size = r1, mu = mean)  # Neg Binomial r = 1.78
#pmf_nb2 <- dnbinom(x, size = r2, mu = mean)  # Neg Binomial r ≈ 0.562
pmf_poisson <- dpois(x, lambda = mean)       # Poisson λ = 0.58

# Create the plot
plot(x, pmf_nb1, type = "h", col = "blue", lwd = 2, 
     ylim = c(0, max(pmf_nb1, pmf_nb2, pmf_poisson) * 1.1),
     xlab = "Number of Events", ylab = "Probability",
     main = paste("Distribution Comparison (mean =", mean, ")"))
points(x, pmf_nb1, col = "blue", pch = 16)

# lines(x, pmf_nb2, type = "h", col = "red", lwd = 2)
# points(x, pmf_nb2, col = "red", pch = 16)

lines(x, pmf_poisson, type = "h", col = "green", lwd = 2)
points(x, pmf_poisson, col = "green", pch = 16)

# Add legend
legend("topright", 
       legend = c(paste("Neg Bin (r =", round(r1, 2), ")"), 
                  #paste("Neg Bin (r =", round(r2, 3), ")"), 
                  paste("Poisson (λ =", mean, ")")),
       col = c("blue",  "green"), lty = 1, lwd = 2)

# Add grid
grid()

