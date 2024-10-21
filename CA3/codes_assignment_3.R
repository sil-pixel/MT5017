
# Call functions built in the assignment 1
source("../Codes/Functions.R")


#########################  ASSIGNMENT_3  #########################
# Load Data
load("../proj_data.Rdata")
modell <- glm(Resultat ~ Alder + Kon + Utbildare, 
              data = data_individ, 
              family = "binomial")
summary(modell)

# Assign values for X and y
y <- matrix(data_individ$Resultat, ncol = 1)
X <- model.matrix(Resultat ~ Alder + Kon + Utbildare, 
                  data = data_individ)

# ---- Task_1 ----

# Compute AIC = 2k - 2l(theta_ml)
n <- nrow(X)
theta0 = rep(0, ncol(X))
k <- length(theta0)
theta_estimate <- NR(theta0, 3, y, X) 
log_likelihood <- l(theta_estimate, y, X)
aic_computed <- 2*k - 2*log_likelihood

#AIC output from R summary 
r_summary_aic <- summary(modell)$aic

# Compute k_cv =sum(l_i(theta_i))/n
# Here, theta_i = theta_ml(X_-i)
nk_cv <- 0
for(i in 1:n) {
  X_minus_i <- X[-i, , drop=FALSE]
  y_minus_i <- y[-i, , drop=FALSE]
  X_i <- X[i, , drop=FALSE]
  y_i <- y[i, , drop=FALSE]
  theta_i <- NR(theta0, 3, y_minus_i, X_minus_i)
  # log likelihood for i-th observation
  log_likelihood_theta_i <- l(theta_i, y_i, X_i)
  nk_cv <- nk_cv + log_likelihood_theta_i
}

k_cv <- nk_cv/n

# Creating a comparison data frame 
comparison_aic_values <- data.frame(
  "AIC_R_model" = r_summary_aic,
  "AIC_computed" = aic_computed,
  "2*nK_CV_computed" = 2*nk_cv
)

# ---- Task_2 ----


# ---- Task_3 ----

