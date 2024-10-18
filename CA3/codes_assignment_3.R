
# Call functions built in the assignment 1
source("../Codes/Functions.R")


#########################  ASSIGNMENT_3  #########################
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


# ---- Task_2 ----


# ---- Task_3 ----

