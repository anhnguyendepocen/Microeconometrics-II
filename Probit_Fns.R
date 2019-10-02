# Shai Knight-Winnig 2017
# Helper Functions: Probit Monte Carlo Simulation

#Define some support functions:
PHI <- function(value){
  pnorm(value, mean = 0, sd = 1) #return CDF of standard normal evaluated at 'value'
}

phi <- function(value){
  dnorm(value, mean = 0, sd = 1) #return PDF of standard normal evaluated at 'value'
}

#Log Likelihood Fn:
Probit_LL <- function(y, x, LL_par) {
  f <- sum(y*log(PHI(x %*% LL_par))) + sum((1-y)*log(1-PHI(x %*% LL_par)))
  return(-f) #Negate here since optim() fn solves min problem by default
}

#Compute probit first order condition for a single row of data (will use for score calc)
#FOC_par is parameter (e.g., beta, gamma, etc.) that you ultimately want to estimate
Probit_FOC <- function(y, x, FOC_par){
  #FOC for probit model (individual row of data):
  ((y - PHI(x %*% FOC_par)) / (PHI(x %*% FOC_par)*(1-PHI(x %*% FOC_par))) * phi(x %*% FOC_par)) %*% x
}

#Gradient/Score (First Order Condition) of Log Likelihood Fn:
#score_par is parameter (e.g., beta, gamma, etc.) that you ultimately want to estimate
Probit_LL_score <- function(y, x, score_par) {
  score <- t(c(0,0)) #initialize score vector (treat it as row for now)
  
  #Loop through each observation, add FOC calc to total sum (this is summing FOC over all i)
  for (i in 1:n){
    score <- score + Probit_FOC(y = y[i], x = x_obs[i,], FOC_par = score_par)
    #cat(i,"Score is:", score, "\n")
  }
  return(-score) #Negate here since optim() fn solves min problem by default 
}


#Calculate Generalized Residual for one observation:
Probit_genResidual_IND <- function(y, x, gen_par){
  ((y - PHI(x %*% gen_par)) / (PHI(x %*% gen_par)*(1-PHI(x %*% gen_par))) * phi(x %*% gen_par))
}

#Calculate average of Generalized Residual across all observations (for one trial):
Probit_genResidual_AVG <- function(y, x, par){
  #For each trial, calculate genResidual by computing for each observation and summing over i
  genResidual <- 0 #initialize generalized residual to 0
  for (i in 1:n){
    genResidual <-  genResidual + Probit_genResidual_IND(y = y[i], x = x[i,], gen_par = par)
      #((y[i] - PHI(x[i,] %*% par)) / (PHI(x[i,] %*% par)*(1-PHI(x[i,] %*% par))) * phi(x[i,] %*% par))
  }
  genResidual <- genResidual / n #divide sum by n to get average genResidual
  return(genResidual)
}

#################################### Tobit Fns ####################################

#Log Likelihood Fn:
Tobit_LL <- function(y, x, D, LL_sigma, LL_par) {
  f <- sum(D*log((phi((y - x %*% LL_par)/LL_sigma))/LL_sigma)) + sum((1-D)*log(1-PHI(x %*% LL_par/LL_sigma)))
  return(-f) #Negate here since optim() fn solves min problem by default
}
