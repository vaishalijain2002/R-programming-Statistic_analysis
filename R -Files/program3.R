#Question 1
microRNA_size <- 20
purine_probability <- 0.7
sample_size <- 100

mu <- microRNA_size * purine_probability
sigma <- sqrt(microRNA_size * purine_probability * (1 - purine_probability))

mu_Y <- mu
sigma_Y <- sigma / sqrt(sample_size)

z <- (15 - mu_Y) / sigma_Y

probability_Y_gt_15 <- 1 - pnorm(z)

cat("Probability that Y (average purines) is greater than 15:", probability_Y_gt_15, "\n")

#Question 2
# Set parameters
X_mean <- 7
X_var <- 3
Y_mean <- 12
Y_variance <- 7
covar <- 3
sample <- 100
simul <- 10000

# Set up variables
exceed_count <- 0
results <- numeric(simul)

# Make simulations
for (i in 1:simul) {
  # Generate samples
  X_samples <- rnorm(sample, mean = X_mean, sd = sqrt(X_var))
  Y_samples <- rnorm(sample, mean = Y_mean, sd = sqrt(Y_variance))
  
  # Determine sample means
  mean_X_sample <- mean(X_samples)
  mean_Y_sample <- mean(Y_samples)
  
  # See if the mean difference is more than 0.5.
  if (mean_Y_sample - mean_X_sample > 0.5) {
    exceed_count <- exceed_count + 1
  }
  results[i] <- mean_Y_sample - mean_X_sample
}

probability <- exceed_count / simul
margin_of_error <- 1.96 * sqrt((probability * (1 - probability)) / simul)
lower_bound <- probability - margin_of_error
upper_bound <- probability + margin_of_error

cat("Probability:", probability, "\n")
cat("95% Confidence Interval: [", lower_bound, ",", upper_bound, "]\n")

#Question 3
num_simulations <- 10000
df_chi_sq <- 8
df_t_dist <- 5

mean_Y_result <- 0

for (iteration in 1:num_simulations) {
  X1 <- sqrt(rchisq(1, df_chi_sq))
  X2 <- rt(1, df_t_dist)
  X3 <- rchisq(1, df_chi_sq)
  
  Y <- X1 * X2 + X3
  
  mean_Y_result <- mean_Y_result + Y
}

mean_Y <- mean_Y_result / num_simulations

cat("Mean of Y:", mean_Y, "\n")

#Question 4
set.seed(123)

sample_size <- 1000
a_n <- sqrt(2 * log(sample_size)) - 0.5 * (log(log(sample_size)) + log(4 * pi)) * (2 * log(sample_size))^(-1/2)
b_n <- (2 * log(sample_size))^(-1/2)

maxima_values <- replicate(1000, max(rnorm(sample_size)))
normalized_maxima <- (maxima_values - a_n) / b_n

extreme_value_function <- function(x) exp(-x) * exp(-exp(-x))
x_values <- seq(0, 5, length.out = 100)

hist(normalized_maxima, prob = TRUE, col = "lightblue", main = "Extreme Value Distribution")
lines(x_values, extreme_value_function(x_values), col = "red", lwd = 2)
curve(dnorm(x), add = TRUE, col = "green", lwd = 2)

legend("topright", legend = c("Normalized Maxima", "Extreme Value Function", "Normal Distribution"),
       col = c("lightblue", "red", "green"), lwd = 2)