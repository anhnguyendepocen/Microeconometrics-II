# Shai Knight-Winnig 2017
#HW 1 - Micro Metrics
#Probit Monte Carlo simulation

rm(list = ls()) #clear environment
library(numDeriv) #needed for grad() function

#Initialize some values/data:
n <- 200 #Number of observations in 'sample'
t <- 100 #Number of trials in simulation
#set.seed(1) #use this if you want to have a "known" sequence of random numbers

#Create x values (which will remain constant througout simulation). One regressor plus an intercept.
x_obs <- cbind(matrix(1, n, 1), matrix(rnorm(n, 0, 1)))

#Generate nXt matrix of random errors (so each trial has slightly different errors):
u_mtx <- matrix(data = rnorm(n*t, mean = 0, sd = 1), nrow = n, ncol = t) 
#u <- rnorm(n, 0, 0.5)

trueBeta <- c(0.2, 0.3) #arbitrarily pick true Beta (in data generating process)

y_star_mtx <- matrix(999, n, t) #create matrix of y_star data (latent model)
y_obs_mtx <- matrix(999, n, t) #create matrix of y_obs (observed y values, will be 0 or 1)

#Populate latent y_star and observed y_obs column by column (trial by trial) with random data
#Use latent model and y_star>0 threshold to populate y_star_mtx and y_obs_mtx 
for (j in 1:t){
  y_star_mtx[,j] <- x_obs %*% trueBeta + u_mtx[,j] #latent model uses x and u to determine y_star
  y_obs_mtx[,j] <- ceiling(y_star_mtx[,j] / max(abs(y_star_mtx))) #ceiling rounds to highest int: neg->0 and pos->1
}

#Prior method used to compute observed y values:
#y_star <- x_obs %*% trueBeta + u
#y_obs <- ceiling(y_star / max(abs(y_star))) #ceiling rounds to highest int: neg->0 and pos->1

#Append observed y values to x values
#dat <- data.frame(x_obs , y_obs) #stores (x,y) data for regression


#Define some support functions:
PHI <- function(value){
  pnorm(value, mean = 0, sd = 1) #return CDF of standard normal evaluated at 'value'
}

phi <- function(value){
  dnorm(value, mean = 0, sd = 1) #return PDF of standard normal evaluated at 'value'
}

#Log Likelihood Fn:
Probit_LL <- function(y, x, beta) {
  f <- sum(y*log(PHI(x %*% beta))) + sum((1-y)*log(1-PHI(x %*% beta)))
  return(-f) #Negate here since optim() fn solves min problem by default
}

#Compute probit first order condition for a single row of data (will use for score calc)
Probit_FOC <- function(y, x, beta){
  #FOC for probit model (individual row of data):
  ((y - PHI(x %*% beta)) / (PHI(x %*% beta)*(1-PHI(x %*% beta))) * phi(x %*% beta)) %*% x
}

#Gradient/Score (First Order Condition) of Log Likelihood Fn:
Probit_LL_score <- function(y, x, beta) {
  score <- t(c(0,0)) #initialize score vector (treat it as row for now)
 
  #Loop through each observation, add FOC calc to total sum (this is summing FOC over all i)
  for (i in 1:n){
    #score <- score + Probit_FOC(y = y_obs[i], x = x_obs[i,], beta = trueBeta)
    score <- score + Probit_FOC(y = y[i], x = x_obs[i,], beta = trueBeta)
    #cat(i,"Score is:", score, "\n")
  }
  return(-score) #Negate here since optim() fn solves min problem by default 
}

#Marginal Effect (Probit):
Probit_ME <- function(x, beta) {
  phi(x %*% beta) %*% beta[2] #'beta[2]' is second element in beta vector (in this case our x coefficient)
}

#Calculate Generalized Residual (for one trial)
Probit_genResidual <- function(y, x, beta){
  #For each trial, calculate genResidual by computing for each observation and summing over i
  genResidual <- 0 #initialize generalized residual to 0
  for (i in 1:n){
    genResidual <-  genResidual + 
      ((y[i] - PHI(x[i,] %*% beta)) / (PHI(x[i,] %*% beta)*(1-PHI(x[i,] %*% beta))) * phi(x[i,] %*% beta))
      #(((y[i] - PHI(x[i,] %*% beta)) / (PHI(x[i,] %*% beta)*(1-PHI(x[i,] %*% beta))) * phi(x[i,] %*% beta)) %*% x[i,])
    }
  return(genResidual)
}

#Check to see if Gradient function is accurate (vs. numerical approximation)
for (z in 1:t){
  grad_check <- grad(function(q) Probit_LL(y = y_obs_mtx[,z] ,x = x_obs, q), trueBeta)
  score_check <- Probit_LL_score(y_obs_mtx[,z], x_obs, trueBeta)
  #print output to verify:
  cat(z, "Score Function Validation:\nScore function output is:", score_check, "\nGradient function output is:", grad_check, "\n")
}

#Solve probit for betaHat using optim() fn.  Note: reltol is convergence tolerance 
betaHat_mtx <- matrix(0, t, 2) #create empty martix to hold betaHat estimates

for (k in 1:t){
  optimization <- optim(par = trueBeta, fn = Probit_LL, y = y_obs_mtx[,k], x = x_obs, gr = Probit_LL_score,
                        method = c("BFGS"), control = list(reltol = 1e-9))

  betaHat_mtx[k,] <- optimization$par #assign correct part of optim() output to betaHat_mtx for kth trial
  cat("betaHat for Trial", k, "is", betaHat_mtx[k,], "\n") #print output to verify
}

#Calculate average values post-simulation
Avg_betaHat <- colMeans(betaHat_mtx)

#Compute Generalized Residuals for all trials and compute average for simulation:
genResidual_mtx <- matrix(999, nrow = 1, ncol = t) #initialize values to 999

#Loop through each trial and compute/store genResidual:
for (h in 1:t){
  genResidual_mtx[1,h] <- Probit_genResidual(y = y_obs_mtx[,h], x = x_obs, beta = betaHat_mtx[h,])
}

#Average across all trials and print output:
Avg_genResidual <- rowMeans(genResidual_mtx)/n

#Calculate Marginal Effects:
#initialize min/max MEs that we'll report at the end: 
minIntercept_ME <- 1e9
maxIntercept_ME <- -1e9
minRegressor_ME <- 1e9
maxRegressor_ME <- -1e9

for (w in 1:t){
  #calculate MEs for each trial and update min/max values
  marginalEffects <- phi(x_obs %*% betaHat_mtx[w,]) %*% t(betaHat_mtx[w,])
  minIntercept_ME <- min(minIntercept_ME, marginalEffects[,1])
  maxIntercept_ME <- max(maxIntercept_ME, marginalEffects[,1])
  minRegressor_ME <- min(minRegressor_ME, marginalEffects[,2])
  maxRegressor_ME <- max(maxRegressor_ME, marginalEffects[,2])
}

#print results from simulations;
cat("Average values for betaHat after", t, "trials are:\n", Avg_betaHat, "\n")
cat("Average value of Generalized Residual after", t, "trials is:", Avg_genResidual, "\n")
cat("Across", t, "trials, the marginal effect for the intercept ranges from", minIntercept_ME, "to", maxIntercept_ME, "\n")
cat("Across", t, "trials, the marginal effect for the x regressor ranges from", minRegressor_ME, "to", maxRegressor_ME, "\n")

