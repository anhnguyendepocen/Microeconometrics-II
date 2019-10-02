# Shai Knight-Winnig 2017
# HW 3 - Micro Metrics
# Tobit Regression with Censored Above Monte Carlo simulation

rm(list = ls()) #clear environment 

#Import functions for Probit operations:
source('Probit_Fns.R') 

#Initialize variables and generate data to be used in all trials:
n <- 200 #Number of observations in 'sample'
t <- 200 #Number of trials in simulation
set.seed(21)  #Uncomment this line if you want to use the same seeded random sequence as someone else to compare results

true_beta <- c(-.1, .3, .4) #arbitrarily pick true beta (data generating process)
sigmasq <- 1 #Set the std dev of the error term, do not change

#Create x data (2 regressors plus an intercept)
x_obs <- cbind(matrix(1, nrow=n, ncol=1), matrix(data=rnorm(n*(length(true_beta)-1), mean=0, sd=2), 
                                                 nrow=n, ncol=length(true_beta)-1))
#Tobit Log-Liklihood Function
tobit_LL <- function(y, x, D, par) {
  
  f <- sum((D%*%phi(y-x%*%par))+(1-D)%*%(1-PHI(x%*%par)))
  
  return(-f) #Negate here since optim() fn solves min problem by default
}
#Create matrices/variables to score key data from each simulation:
naive_beta_hat_mtx <- matrix(0, nrow=length(true_beta), ncol=t) #store initial OLS coefficients (used to seed optim fn)
ols_beta_hat_mtx <- <- matrix(0, nrow=length(true_beta), ncol=t)
beta_hat_mtx <- matrix(0, nrow=length(true_beta), ncol=t)

############## MAIN LOOP FOR TOBIT ####################

for (i in 1:t){
  
  ########################## STEP 1: GENERATE DATA FOR SIMULATION ##########################
  
  #Generate the error term
  u <- matrix(data=rnorm(n, mean=0, sd=sigmasq), nrow=n, ncol=1)
  
  #Generate y_star (using true model) and y_obs (using indicator function)
  y_star <- x_obs %*% true_beta + u #latent process to determine y_star
  y_obs <- matrix(0, nrow=n, ncol=1)
  data_cens <- matrix(0, nrow=1, ncol=n)
  for (j in 1:n){
    y_obs[j] <- ifelse(y_star[j]>0,y_star[j], 0) #indicator determines whether y_star is observed 
    data_cens[j]<- ifelse(y_obs[j]==0,0,1)
  }
  
  ###########################################################################################
  
  ###################### STEP 2: RUN TOBIT REGRESSION (MLE) #################################
  
  ###########################################################################################
  
  #Run naive OLS to get seed for optimization function
  naive_beta_hat_mtx[,i] <- solve(t(x_obs) %*% x_obs) %*% t(x_obs) %*% y_obs
  
  #Do Tobit to get estimate of beta
  optimization <- optim(par = naive_beta_hat_mtx[,i], fn = tobit_LL, y = y_obs, x = x_obs, D = data_cens, 
                        method = "Nelder-Mead", hessian = FALSE, control = list(reltol = 1e-9))
  
  beta_hat_mtx[,i] <- optimization$par
  
  
  
  ###########################################################################################
  
  ###################### STEP 3: RUN OLS REGRESSION WITH GEN. RESIDUAL #####################
  
  ###########################################################################################
  
  
  
  
  
  
}

#Output Checks
cat('The number of censored data points are', sum(data_cens),'out of', n, '\n')
cat('The avg beta hat estimate is', rowMeans(beta_hat_mtx), '\n')