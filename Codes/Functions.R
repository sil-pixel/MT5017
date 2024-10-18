
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
I <- function(theta, X){
  p <- 1/(1 + exp(-X %*% theta))
  v <- p*(1-p)
  D <- diag(as.vector(v))
  info <- t(X) %*% D %*% X     # dim: k*k
  return(info)
} 


#---- Task-1 b) Use Newton-Raphson's algorithm to compute the ML-estimates in a logistic regression model ----

NR <- function(theta0, niter, y, X) {
  # Initialize theta with the starting value
  theta <- theta0  
  
  for (i in 1:niter) {
    score <- S(theta, y, X)  
    info <- I(theta, X)   
    
    # Update theta using the Newton-Raphson formula
    theta <- theta + solve(info) %*% score
    
    # Print current estimate (for tracking)
    #cat("Iteration", i, ": theta =", theta, "\n")
  }
  
  return(theta)  # Return the final estimate
}
