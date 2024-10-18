# ---- Task-1 a) Likelihood, log-likelihood, score function, and Fisher information ----
# Dimensionality:
# X: N*k
# p: N*1
# y: N*1
# theta: k*1
# D: N*N
# 
# %*% is the dot product operation


# Likelihood
L <- function(theta, y, X){
  # Compute the linear predictor
  eta <- X %*% theta   # dim: N*1
  
  # Compute the probabilities using the logistic function
  p <- 1/(1 + exp(-eta))   # dim: N*1
  
  likelihood <- prod(dbinom(y, size=1, prob=p))  # scalar
  return(likelihood)
}


# Log-likelihood
l <- function(theta, y, X){
  likelihood <- L(theta, y, X)
  log_likelihood <- log(likelihood)
  return(log_likelihood)
}


# Score function
S <- function(theta, y, X){
  p <- 1/(1 + exp(-X %*% theta))
  score <- t(X) %*% (y-p)   # dim: k*1
  return(as.vector(score))
}


# Fisher information
I <- function(theta, y, X){
  p <- 1/(1 + exp(-X %*% theta))
  v <- p*(1-p)
  D <- diag(as.vector(v))
  info <- t(X) %*% D %*% X     # dim: k*k
  return(info)
} 


#---- Task-1 b) Use Newton-Raphson's algorithm to compute the ML-estimates 
# in a logistic regression model ----

NR <- function(theta0, niter, y, X) {
  # Initialize theta with the starting value
  theta <- theta0  

  for (i in 1:niter) {
    score <- S(theta, y, X)  
    info <- I(theta, y, X)   

    # Update theta using the Newton-Raphson formula
    theta <- theta + solve(info) %*% score

    # Print current estimate (for tracking)
    cat("Iteration", i, ": theta =", theta, "\n")
  }

  return(theta)  # Return the final estimate
}



#---- Load data for task 2 ----
library(RCurl)
set.seed(900101)
länk <- "https://raw.githubusercontent.com/mskoldSU/MT5003_HT17/master/Projekt/proj_data.csv"
data_individ <- read.csv(text = getURL(länk))
idx <- sample(1:nrow(data_individ), 1000)
data_individ <- data_individ[idx, ]
save(data_individ, file = "proj_data.Rdata")


#---- Task-2 ----
modell <- glm(Resultat ~ Alder + Kon + Utbildare, data = data_individ, family = "binomial")
summary(modell)

y <- matrix(data_individ$Resultat, ncol=1)
X <- model.matrix(Resultat ~ Alder + Kon + Utbildare, data = data_individ)
theta0 = c(0, 0, 0, 0)

theta4 <- NR(theta0, 4, y, X)
compare_theta4 <- data.frame("theta_NR4" = theta4, "theta_R_model" = coef(modell))

theta3 <- NR(theta0, 3, y, X) # We will use this value in the following steps.
compare_theta3 <- data.frame("theta_NR3" = theta3, "theta_R_model" = coef(modell))

theta2 <- NR(theta0, 2, y, X)
compare_theta2 <- data.frame("theta_NR2" = theta2, "theta_R_model" = coef(modell))


#---- Task-3 ----
# Compare the standard error values with R 

info <- I(theta3, y, X) #dim = k*k
info_inv <- solve(info)   #dim = k*k
diagonal <- diag(info_inv) # dim = k*1
standard_errors <- sqrt(diagonal) #dim = k*1


# Standard errors from R in summary 
r_se <- summary(modell)$coefficients[, "Std. Error"]


# Creating a comparison data frame 
comparison_se <- data.frame(
  "Se_R_model" = r_se,
  "Se_computed" = standard_errors
)

#---- Task-4 ----

# Bootstrap function for the SE
bootstrappingSE <- function(f_hat, n){
  # Matrix to store the bootstrap estimates for each coefficient
  bootstrap_estimates <- matrix(NA, nrow = n, ncol = length(coef(f_hat)))
  
  for (i in 1:n) {
    # Step 1: Get the predicted probabilities
    p_hat <- predict(f_hat, type = "response")
    
    # Step 2: Resample y based on the predicted probabilities
    y_boot <- rbinom(n = length(p_hat), size = 1, prob = p_hat)
    
    # Step 3: Refit the logistic model using the bootstrapped y
    model_boot <- glm(y_boot ~ Alder + Kon + Utbildare, data = data_individ, family = binomial)
    
    # Step 4: Store the bootstrapped estimates
    bootstrap_estimates[i, ] <- coef(model_boot)
  }
  return(bootstrap_estimates)
}

# Calculate Standard Errors
n_bootstrap <- 10000
bootstrap_estimates <- bootstrappingSE(modell, n_bootstrap)

# Calculate the standard deviation of the bootstrapped coefficients (standard errors)
bootstrap_se <- apply(bootstrap_estimates, 2, sd)

# Add the bootstrap standard errors to the comparison matrix
comparison_se$Se_bootstrap <- bootstrap_se


# Bootstrap 95% confidence interval for the probability that someone privately educated,
#of age 34, female is successful
bootstrappingP <- function(f_hat, new_data, n){
  # Vector to store the bootstrap probabilities
  p_bootstrap <- numeric(n)
  
  for (i in 1:n) {
    # Step 1: Get the predicted probabilities
    p_hat <- predict(f_hat, type = "response")
    
    # Step 2: Resample y based on the predicted probabilities
    y_boot <- rbinom(n = length(p_hat), size = 1, prob = p_hat)
    
    # Step 3: Refit the logistic model using the bootstrapped y
    model_boot <- glm(y_boot ~ Alder + Kon + Utbildare, data = data_individ, family = binomial)
    
    # Step 4: Predict for the specific individual using the bootstrapped model
    p_bootstrap[i] <- predict(model_boot, newdata = new_data, type = "response")
  }
  return(p_bootstrap)
}


X_new <- data.frame("Alder" = 34, "Kon"="Kvinna", "Utbildare"="Privatist")
n_bootstrap_p <- 10000
p_bootstrap <- bootstrappingP(modell, X_new, n_bootstrap_p)

# Calculate the 95% Confidence Interval
bootstrap_ci <- quantile(p_bootstrap, c(0.025, 0.975))
