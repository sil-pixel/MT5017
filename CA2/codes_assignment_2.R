
# Call functions built in the assignment 1
source("Functions.R")


#########################  ASSIGNMENT_2  #########################
# Load Data
load("proj_data.Rdata")
modell <- glm(Resultat ~ Alder + Kon + Utbildare, 
              data = data_individ, 
              family = "binomial")
summary(modell)

# Assign values for X and y
y <- matrix(data_individ$Resultat, ncol = 1)
X <- model.matrix(Resultat ~ Alder + Kon + Utbildare, 
                  data = data_individ)

# ---- Task_1 ----

#Verify if z values are Wald statistics using I and NR functions

# Wald statistic = (theta_estimate - theta0)/(standard_error)
#Consider the theta estimate after 3 iterations as per Part-I
theta0 = c(0, 0, 0, 0)
theta_estimate <- NR(theta0, 3, y, X) 
I_matrix <- I(theta_estimate, X)
standard_error <- sqrt(diag(solve(I_matrix)))
wald_statistic <- theta_estimate/standard_error

#Compare the wald estimate output with z values in summary 
z_values <- summary(modell)$coefficients[, "z value"]

# Creating a comparison data frame 
comparison_z_values <- data.frame(
  "Z_values_R_model" = z_values,
  "Z_values_computed" = wald_statistic
  )

# ---- Task_2 ----

# Compute generalized likelihood ratio statistics that corresponds 
# to Wald Statistics in Task1
# GN_L = 2 x log(profile_likelihood(theta_estimate)/profile_likelihood(theta0))

# Generalized likelihood ratio statistics
compute_gnl_ratio <- function(XP) {
  #compute profile log likelihood for theta estimate
  profile_theta <- rep(0, ncol(XP))
  profile_theta_estimate <- NR(profile_theta, 3, y, XP)
  profile_log_likelihood_estimate <- l(profile_theta_estimate, y, XP)
  
  theta_estimate0 <- NR(theta0, 3, y, X)  #### !!!! different length (4) than profile_theta_estimate (3)!!!
  log_likelihood_estimate0 <- l(theta_estimate0, y, X)
  profile_theta0 <- rep(0, ncol(XP))
  profile_theta_estimate0 <- NR(profile_theta0, 3, y, XP)
  profile_log_likelihood_estimate0 <- l(profile_theta_estimate0, y, XP)
  
  gnl_ratio <- 2* (log_likelihood_estimate0 - profile_log_likelihood_estimate)
  # theta_estimate from task1 generated for Wald Statistics since it
  # contains the MLEs for theta and eta
  log_likelihood_estimate <- l(theta_estimate, y, X)

  gnl_ratio <- 2* (log_likelihood_estimate - profile_log_likelihood_estimate0)
  return(gnl_ratio)
}

#compute the gnl ratio for each parameter and profile the X matrix by each parameter
gnl_ratio_intercept <- compute_gnl_ratio(X[, -1])
gnl_ratio_alder <- compute_gnl_ratio(X[, -2])
gnl_ratio_kon <- compute_gnl_ratio(X[ ,-3])
gnl_ratio_utbildare <- compute_gnl_ratio(X[, -4])

# Creating a comparison data frame for gnl ratio
comparison_gnl_values <- data.frame(
  "squared_wald_statistic" = wald_statistic^2,
  "generalized_likelihood_ratio_statistics" = 
    c(gnl_ratio_intercept, gnl_ratio_alder, gnl_ratio_kon, gnl_ratio_utbildare)
)

# Determine the corresponding P-values
p_value_intercept <- pchisq(gnl_ratio_intercept, df = 1, lower.tail = FALSE)
p_value_alder <- pchisq(gnl_ratio_alder, df = 1, lower.tail = FALSE)
p_value_kon <- pchisq(gnl_ratio_kon, df = 1, lower.tail = FALSE)
p_value_utbildare <- pchisq(gnl_ratio_utbildare, df = 1, lower.tail = FALSE)

p_values <- data.frame(
  "intercept_profile" = p_value_intercept,
  "alder_profile" = p_value_alder, 
  "kon_profile" = p_value_kon, 
  "utbildare_profile" = p_value_utbildare)

#Compare the p values output with p values in summary 
p_values_summary <- summary(modell)$coefficients[, "Pr(>|z|)"]

# Creating a comparison data frame 
comparison_p_values <- data.frame(
  "P_values_R_model" = p_values_summary,
  "P_values_computed" = c(p_value_intercept, p_value_alder, p_value_kon, p_value_utbildare)
)



# ---- Task_3 ----

# Take columns "Alder" and "Utbildare", drop the rest (assuming parameter of 0 due to H_0: theta = theta0 = (0, 0))
X_new = X[, c(2,4) ]

# Estimate the parameter vector eta for the two columns
eta <- NR(theta0=c(0, 0), niter=10, y, X_new) 

# Given eta and theta0, compute the generalized score statistic Ts_theta0
theta_new <- c(0, eta[1], 0, eta[2])
score <- S(theta_new, y, X)
info <- I(theta_new, X)
info_inverse <- solve(info)
Ts_theta0 <- t(score) %*% info_inverse %*% score

# Calculate the p-value
df <- length(theta_new)  # Degree of freedom; number of parameters under H_0
p_value <- 1 - pchisq(Ts_theta0, df = df)


# ---- Task_4 ----

# Grid values for theta_Kon
grid <- seq(0.4-1, 0.4+1, 0.01)

# Compute the profile likelihood
L_p <- numeric(length(grid))
L_e <- numeric(length(grid)) 

for(i in seq(length(grid))){
  theta_Kon <- grid[i] 
  new_model <- glm.fit(x = X[, -3], y = y,
                  offset = theta_Kon * X[, 3],
                  family = binomial())
  theta_new <- c(new_model$coeff[1:2], theta_Kon, new_model$coeff[3])
  
  # Compute the profile likelihood
  L_p[i] <- L(theta_new, y, X)
  
  # Compute the estimated likelihood
  L_e[i] <- L(c(theta_estimate[1:2], theta_Kon, theta_estimate[4]), y, X)
}


# Plot the profile likelihood and the estimated likelihood
## Create the plot
plot(grid, L_p, type = "l", col = "darkgreen", lwd=3,
     main = "Profile and Estimated Likelihood of theta_Kon",
     xlab = "theta_Kon",
     ylab = "Likelihood")

## Add the estimated likelihood
lines(grid, L_e, col = "blue", lty = 4, lwd=3)

## Draw a horizontal line for the 95% profile likelihood confidence interval
ci_95 <- max(L_p) * exp(-0.5 * qchisq(0.95, df=1))
abline(h = ci_95, col = "red", lty = 3, lwd=3) 

## Add a legend
legend("topright", legend = c("Profile likelihood", "Estimated likelihood", "95% threshold"),
       col = c("darkgreen", "blue", "red"), lty = c(1, 4, 3), lwd = 3)


# Compare with the Wald confidence interval
se_Kon <- standard_error[3]
theta_Kon_estimated <- theta_estimate[3]
lb <- theta_Kon_estimated - qnorm(0.975) * se_Kon
ub <- theta_Kon_estimated + qnorm(0.975) * se_Kon
cat(sprintf("[%f, %f]", lb, ub))

