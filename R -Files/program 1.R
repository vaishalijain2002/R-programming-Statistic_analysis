#problem 2 
n <- 1000
result <- n*((n + 1) * (2 * n + 1)) / 6
print(paste("The result of ∑(x=1)^1000 x^2 =", result))

#problem 3(a)
age <- c(21, 25, 31, 35, 41, 45, 51, 55, 61, 65)
#problem 3(b)
age_converts_to_months <- age * 12
#problem 3(c)
total_sum_of_age_of_all_months <- sum(age_converts_to_months)
#problem 3 (d)
youngest_among_of_all <- min(age)
#problem 3 (e)
oldest_among_of_all <- max(age)
#problem 3 (f)
squareroot_of_age <- sqrt(age)


#problem 4(g)
X <- 3 * (1:30)
print(X)

#problem 4(h)
Y <- rep(0, 30)
cat(" Y:", Y, "\n")

#problem 4(i)
the_y_vector <- rep(0, 30)
for (k in 1:30) {
  the_y_vector[k] <- if (k < 20) sin(2 * k) else integrate(function(x) x, lower = 0, upper = k)$value
}
cat("Vector Y:", the_y_vector, "\n")




