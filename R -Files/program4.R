#Question1
#a) Observed data
dist_results <- c(1.433, 0.524, 0.384, 4.515, 1.852, 0.429)

# Likelihood function
likelihood <- function(lambda) {
  -sum(log(dexp(dist_results, rate = lambda)))
}

# Initial bounds for lambda
lower_bound <- 0.001
upper_bound <- 10

# Optimization using "optimize"
result <- optimize(f = likelihood, interval = c(lower_bound, upper_bound))

# MLE estimate
mle_numerical <- result$minimum
mle_numerical

#b)# Data
data <- c(1.433, 0.524, 0.384, 4.515, 1.852, 0.429)

# Calculate MLE for lambda
mle_lambda <- 6 / sum(data)

# Print the MLE
print(paste("MLE for lambda:", mle_lambda))


#Question 2 
# Function to calculate the method of moments estimator of m
calculate_moments_estimator <- function(sample_mean) {
  # The method of moments estimator for m is the sample mean itself
  estimator_m <- sample_mean
  return(estimator_m)
}

# Function to calculate the lower confidence interval for m
calculate_lower_confidence_interval <- function(sample_mean, sample_sd, n, alpha = 0.1) {
  # Calculate the t-score for a one-sided 90% confidence interval
  t_score <- qt(1 - alpha, df = n - 1)
  
  # Calculate the margin of error
  margin_of_error <- t_score * (sample_sd / sqrt(n))
  
  # Calculate the lower bound of the confidence interval
  lower_bound <- sample_mean - margin_of_error
  
  return(lower_bound)
}

# Given values (replace with your actual data)
sample_mean <- 98.6
sample_sd <- 9.4
n <- 75

# (a) Calculate the method of moments estimator of m
moment_estimator <- calculate_moments_estimator(sample_mean)
cat("Method of Moments Estimator of m:", moment_estimator, "\n")

# (b) Calculate the one-sided 90% Lower Confidence Interval of m
lower_ci <- calculate_lower_confidence_interval(sample_mean, sample_sd, n)
cat("One-sided 90% Lower Confidence Interval of m:", lower_ci, "\n")

#Question 3
# Calculate bootstrap 95% CIs for the mean of gene expression in the "ALL" group
mean_ci_all_bootstrap <- quantile(boot_means_all, c(0.025, 0.975))

# Calculate bootstrap 95% CIs for the variance of gene expression in the "ALL" group
variance_ci_all_bootstrap <- quantile(boot_variances_all, c(0.025, 0.975))

# Calculate bootstrap 95% CIs for the mean of gene expression in the "AML" group
mean_ci_aml_bootstrap <- quantile(boot_means_aml, c(0.025, 0.975))

# Calculate bootstrap 95% CIs for the variance of gene expression in the "AML" group
variance_ci_aml_bootstrap <- quantile(boot_variances_aml, c(0.025, 0.975))

# Print the results
cat("Bootstrap 95% CI for Mean (ALL):", mean_ci_all_bootstrap, "\n")
cat("Bootstrap 95% CI for Variance (ALL):", variance_ci_all_bootstrap, "\n")
cat("Bootstrap 95% CI for Mean (AML):", mean_ci_aml_bootstrap, "\n")
cat("Bootstrap 95% CI for Variance (AML):", variance_ci_aml_bootstrap, "\n")
# Calculate t-intervals for the mean of gene expression in the "ALL" group
t_interval_all <- mean(all_group) + qt(c(0.025, 0.975), df = n - 1) * sd(all_group) / sqrt(n)

# Calculate t-intervals for the mean of gene expression in the "AML" group
t_interval_aml <- mean(aml_group) + qt(c(0.025, 0.975), df = length(aml_group) - 1) * sd(aml_group) / sqrt(length(aml_group))

# Print the results
cat("T-Interval 95% CI for Mean (ALL):", t_interval_all, "\n")
cat("T-Interval 95% CI for Mean (AML):", t_interval_aml, "\n")
# Calculate bootstrap 95% CI for the median gene expression in the "ALL" group
median_ci_all_bootstrap <- quantile(boot_means_all, c(0.025, 0.975))

# Calculate bootstrap 95% CI for the median gene expression in the "AML" group
median_ci_aml_bootstrap <- quantile(boot_means_aml, c(0.025, 0.975))

# Print the results
cat("Bootstrap 95% CI for Median (ALL):", median_ci_all_bootstrap, "\n")
cat("Bootstrap 95% CI for Median (AML):", median_ci_aml_bootstrap, "\n")

#Question 4
# Function to calculate the Poisson confidence interval based on sample mean
poisson_mean_ci <- function(data, alpha = 0.1) {
  n <- length(data)
  lambda_hat <- mean(data)
  margin <- qnorm(1 - alpha/2) * sqrt(lambda_hat / n)
  lower_bound <- lambda_hat - margin
  upper_bound <- lambda_hat + margin
  return(c(lower_bound, upper_bound))
}

# Function to calculate the Poisson confidence interval based on sample variance
poisson_variance_ci <- function(data, alpha = 0.1) {
  n <- length(data)
  lambda_hat <- mean(data)
  var_hat <- var(data)
  chi2 <- qchisq(c(alpha/2, 1 - alpha/2), df = n - 1)
  lower_bound <- (n - 1) * var_hat / chi2[2]
  upper_bound <- (n - 1) * var_hat / chi2[1]
  return(c(lower_bound, upper_bound))
}

monte_carlo_simulation <- function(lambda, nsim = 1000, sample_size = 50) {
  coverage_mean <- coverage_variance <- numeric(nsim)
  
  # Generate nsim datasets from the Poisson distribution
  for (i in 1:nsim) {
    data <- rpois(sample_size, lambda)
    
    # Calculate CIs using the updated function names
    ci_mean <- poisson_mean_ci(data)
    ci_variance <- poisson_variance_ci(data)
    
    # Check if true lambda is within CIs
    coverage_mean[i] <- lambda >= ci_mean[1] && lambda <= ci_mean[2]
    coverage_variance[i] <- lambda >= ci_variance[1] && lambda <= ci_variance[2]
  }
  
  # Return proportion of coverage for both CIs
  return(c(mean(coverage_mean), mean(coverage_variance)))
}

# Run simulation for lambda values 0.1, 1, and 10
lambda_values <- c(0.1, 1, 10)
nsim <- 1000

results <- matrix(0, nrow = length(lambda_values), ncol = 2)

# Run the Monte Carlo simulation for nsim runs at three different parameter values
for (i in 1:length(lambda_values)) {
  results[i,] <- monte_carlo_simulation(lambda_values[i], nsim)
}

colnames(results) <- c("Coverage Mean CI", "Coverage Variance CI")
rownames(results) <- paste("Lambda =", lambda_values)
print(results)


