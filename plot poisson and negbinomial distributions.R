

# Set seed for reproducibility
set.seed(1234)

# Parameters
mean <- 0.58
r1 <- 1.78          # First dispersion parameter

# Generate 1000 random variates
n <- 100000
poisson_rvs <- rpois(n, lambda = mean)
nb1_rvs <- rnbinom(n, size = r1, mu = mean)

# Define x range for plotting (based on max observed values)
x_max <- max(poisson_rvs, nb1_rvs) #, nb2_rvs)
x <- 0:x_max

# Calculate theoretical PMFs
pmf_poisson <- dpois(x, lambda = mean)
pmf_nb1 <- dnbinom(x, size = r1, mu = mean)
#pmf_nb2 <- dnbinom(x, size = r2, mu = mean)

# Create histograms with relative frequencies and labels
par(mfrow = c(2, 1))  # 3 plots in a vertical stack

# Poisson
hist_data_p <- hist(poisson_rvs, breaks = seq(-0.5, x_max + 0.5, 1), freq = FALSE, 
                    col = rgb(0, 1, 0, 0.5), main = paste("Poisson (λ =", mean, ")"), 
                    xlab = "Value", ylab = "Relative Frequency", ylim = c(0, 0.9))
counts_p <- hist(poisson_rvs, breaks = seq(-0.5, x_max + 0.5, 1), plot = FALSE)$counts
percent_p <- round(counts_p / n * 100, 1)
labels_p <- paste(counts_p, "\n(", percent_p, "%)", sep = "")
text(x, hist_data_p$density + 0.05, labels = labels_p, pos = 3, col = "green")
points(x, pmf_poisson, col = "green", pch = 16, cex = 1.5)
lines(x, pmf_poisson, col = "green", lwd = 2)

# Negative Binomial (r = 1.78)
hist_data_nb1 <- hist(nb1_rvs, breaks = seq(-0.5, x_max + 0.5, 1), freq = FALSE, 
                      col = rgb(0, 0, 1, 0.5), main = paste("Neg Bin (r =", round(r1, 2), ")"), 
                      xlab = "Value", ylab = "Relative Frequency", ylim = c(0, 0.9))
counts_nb1 <- hist(nb1_rvs, breaks = seq(-0.5, x_max + 0.5, 1), plot = FALSE)$counts
percent_nb1 <- round(counts_nb1 / n * 100, 1)
labels_nb1 <- paste(counts_nb1, "\n(", percent_nb1, "%)", sep = "")
text(x, hist_data_nb1$density + 0.05, labels = labels_nb1, pos = 3, col = "blue")
points(x, pmf_nb1, col = "blue", pch = 16, cex = 1.5)
lines(x, pmf_nb1, col = "blue", lwd = 2)

# Negative Binomial (r ≈ 0.562)
# hist_data_nb2 <- hist(nb2_rvs, breaks = seq(-0.5, x_max + 0.5, 1), freq = FALSE, 
#                       col = rgb(1, 0, 0, 0.5), main = paste("Neg Bin (r =", round(r2, 3), ")"), 
#                       xlab = "Value", ylab = "Relative Frequency", ylim = c(0, 0.9))
# counts_nb2 <- hist(nb2_rvs, breaks = seq(-0.5, x_max + 0.5, 1), plot = FALSE)$counts
# percent_nb2 <- round(counts_nb2 / n * 100, 1)
# labels_nb2 <- paste(counts_nb2, "\n(", percent_nb2, "%)", sep = "")
# text(x, hist_data_nb2$density + 0.05, labels = labels_nb2, pos = 3, col = "red")
# points(x, pmf_nb2, col = "red", pch = 16, cex = 1.5)
# lines(x, pmf_nb2, col = "red", lwd = 2)

# Reset plotting parameters
par(mfrow = c(1, 1))
