# Shai Knight-Winnig 2017
# Ichimura/Robinson 2 Step Estimation

rm(list = ls()) #clear environment
set.seed(187)

#Importing function for Probit 
source("Fxns_HW4.R") 


#Initialize variables and generate data to be used in all trials:
n <- 200 #Number of observations in 'sample'
t <- 100 #Number of trials in simulation

true_beta <- c(-3, 1, 2) #arbitrarily pick true beta (data generating process)
true_theta <- c(-2, 0.5, 1) #arbitrarily pick true theta (data generating process)


#Create x data to be used in each simulation's last step regression (intercept plus 2 regressors) #OUTCOME EQUATION 
x_obs <- cbind(matrix(data=rnorm(n*(length(true_beta)), mean=0, sd=1), 
                                                 nrow=n, ncol=length(true_beta)))

#Create data to be used in each simulation
z_data <- x_obs 


#Create matrices/variables to score key data from each simulation:
ols_guess_mtx <- matrix(0, nrow=length(true_theta), ncol=t) #store initial OLS coefficients (used to seed optim fn)
ichi_coeff_mtx <- matrix(0, nrow=length(true_theta), ncol=t) #result of step Tobit regression
avg_censored <- 0 #track number of censored data points throughout the simulation.
robinson_OLS <- matrix(0, nrow= length(true_beta), ncol=t)
##Initialize some matrix to store data for nonpar

fxln <- matrix(0,n,1)

########### LOOP FOR ICHIMURA STARTS HERE ###############################################

for (i in 1:t){
  
  ##### GENERATE THE MODEL (OUTCOME EQUATION, SELECTION EQUATION)
  
  #Jointly uniform distributed errors (v is a linear transformation of u)
  v <- matrix(data = runif(n, min=-1, max=1), nrow = n, ncol = 1) 
  u <- 0.5*v + matrix(data = runif(n, min=0, max=1)) ##introduces corr btwn error terms (non-normal)
  
  
  #Generate y_star (wage/outcome equation) 
  y_star <- x_obs %*% true_beta + u #(outcome equation)
  
  #Generate h_star (desired work hours/selection equation) and use 'Indicator Function' via ceiling command:
  h_star <- z_data %*% true_theta + v #(selection equation)
  d_indicator <- ifelse(h_star > 0, 1, 0) #assigns 1 if h*>0, 0 if h*<=0 (indicator fn)
  avg_censored <- avg_censored + n-sum(d_indicator) #keep running sum of censored data   
  
  
  #Generate h_obs (observed when hours>0)
  h_obs <-matrix(999,nrow=n, ncol=1)
  matrix
  for (j in 1:n) {
    
    h_obs[j] <- d_indicator[j]*h_star[j]} #indicator determines whether we osberve h_star or not
  
  
  #Generate y_obs (observed when hours>0)
  y_obs <-matrix(999,nrow=n, ncol=1)
  matrix
  for (j in 1:n) {
    
    y_obs[j] <- d_indicator[j]*y_star[j]} #indicator determines whether we osberve h_star or not
  
  
  
  ##Get rid of censored data
  z_data_cens <- cbind(d_indicator, z_data)
  z_data_cens <- z_data_cens[z_data_cens[,1]==1,] #Keeps all rows where first column equals 1
  z_data_cens <- z_data_cens[,-1] #delete first column (indicator fn no longer needed)
  
  h_obs_cens <- cbind(d_indicator, h_obs)
  h_obs_cens <- h_obs_cens[h_obs_cens[,1]==1,] #Keeps all rows where first column equals 1
  h_obs_cens <- h_obs_cens[,-1] #delete first column (indicator fn no longer needed)
  
  y_data_cens <- cbind(d_indicator, y_obs)
  y_data_cens <- y_data_cens[y_data_cens[,1]==1,] #Keeps all rows where first column equals 1
  y_data_cens <- y_data_cens[,-1] #delete first column (indicator fn no longer needed)
  y_data_cens <- matrix(data=y_data_cens, nrow=sum(d_indicator), ncol=1) #convert into matrix
  
  x_data_cens <- cbind(d_indicator, x_obs)
  x_data_cens <- x_data_cens[x_data_cens[,1]==1,] #Keeps all rows where first column equals 1
  x_data_cens <- x_data_cens[,-1] #delete first column (indicator fn no longer needed)
  
  ##Step 1: Guess beta vector to seed the optim
  ols_guess_mtx[,i] <- solve(t(z_data_cens) %*% z_data_cens) %*% t(z_data_cens) %*% h_obs_cens
  
  #######################################################################################
  ##Step 2: Use this guess, the non-censored z data, and non censored observed h_star to find thetahat
  #######################################################################################
  
  ##This codes the obj funtion: sum((y-E[y|Z*theta])^2)
  objfxn <- function(x, y, par) {
    
    
  z_theta <- x%*%par
  sx <- sd(z_theta) #std dev of our random variable 
  bin_width <- 1.06 #the higher the larger the bin width and the "smoother" the function this value is optimal for Gaussian
  fxn <- matrix(0,length(y),1)
  h <- bin_width*sx*(n^(-.2)) #"bin" size, basically the width of our kernel, this is somehow the optimal value for Gaussian
  for (k in 1:length(fxn)){ #loop should be length equal to "uncensored" sample, not "n".
    
    ##This calcs the "z-score" of each x with 0 at x_i
    d = (z_theta[k]-z_theta)/h
    
    ##Builds the kernel using the std normal pdf and inputing the z score of each x
    fx1 <- dnorm(d, mean = 0, sd = 1)
    num <- fx1*y
    fxn[k] <- colSums(num)/colSums(fx1)
  }
  sumobj <- 0
  for (j in 1:length(fxn)){
    sumobj <- sumobj + (y[j]-fxn[j])^2
  }
  
  f <- (sumobj)

  return(f)
  }
  
  
  ##Step 3: Find thetas that minimize sum{(y-E[y|Z*theta])^2}
  
  #NOTE: Optimization should be run on subsample (uncensored data)
  
  optimization <- optim(par = ols_guess_mtx[,i], fn = objfxn, y = h_obs_cens, x=z_data_cens, 
                        method = "Nelder-Mead", hessian = FALSE)
  
  ichi_coeff_mtx[,i] <- optimization$par
  
  
  ##Step 4: Compute residuals from Ichimura step
  v_hat <- h_obs_cens - z_data_cens %*% ichi_coeff_mtx[,i]
  
  ##Step 5: Compute E[y|v_hat] and E[x|v_hat] for uncensored data (fn defined in 'Probit_Fns.R')
  cond_exp_y <- Nonpar_Cond_Exp(y=y_data_cens, x=v_hat)
  cond_expx1 <- Nonpar_Cond_Exp(y = x_data_cens[,1], x = v_hat)
  cond_expx2 <- Nonpar_Cond_Exp(y = x_data_cens[,2], x = v_hat)
  cond_expx3 <- Nonpar_Cond_Exp(y = x_data_cens[,3], x = v_hat)
  
  cond_exp_x <- cbind(cond_expx1,cond_expx2,cond_expx3)
  
  robinson_OLS <- solve(t(x_data_cens-cond_exp_x)%*%(x_data_cens-cond_exp_x))%*%(t(x_data_cens-cond_exp_x)%*%y_data_cens)
  ##Step 6: Compute 2nd stage OLS for Robinson estimator to get unbiased betas:
  
}

cat(avg_censored/t,"of 200 were censored on average", "\n")
cat("Avg. Theta Hat",rowMeans(ichi_coeff_mtx), "\n")
cat("Avg. Beta Hat",rowMeans(robinson_OLS), "\n")

  
  
  
  
  
  
  
  
  
