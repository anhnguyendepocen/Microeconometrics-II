# Shai Knight-Winnig 2017
# HW 2 - Micro Metrics
# Two Step Sample Selection Monte Carlo simulation

rm(list = ls()) #clear environment 

#Import functions for Probit operations:
source('Probit_Fns.R') 

#Initialize variables and generate data to be used in all trials:
n <- 200 #Number of observations in 'sample'
t <- 500 #Number of trials in simulation
#set.seed(21)  #Uncomment this line if you want to use the same seeded random sequence as someone else to compare results

true_gamma <- c(-0.2, 0.9, -0.7, 0.4, 0.5) #arbitrarily pick true gamma (data generating process)
#true_gamma <- c(-0.2, 0.9, -0.7) #FOR PART 2 OF HW WHEN X=Z (NOTE: true_beta and true_gamma must be same length in this case)
true_beta <- c(-0.6, 0.5, 0.2) #arbitrarily pick true beta (data generating process)

#NOTE: x_obs and z_data will automatically adjust to accommodate the size of the true_gamma and true_beta vectors
#Create x data to be used in each simulation's 2nd step regression (2 regressors plus an intercept)
x_obs <- cbind(matrix(1, nrow=n, ncol=1), matrix(data=rnorm(n*(length(true_beta)-1), mean=0, sd=1), 
                                                 nrow=n, ncol=length(true_beta)-1))

#Create z data to be used in each simulation's 1st step Probit (x_obs plus one more regressor)
z_data <- cbind(x_obs, matrix(data=rnorm(n*(length(true_gamma)-length(true_beta)), mean=0, sd=1), 
                              nrow=n, ncol=length(true_gamma)-length(true_beta)))
#z_data <- x_obs #FOR PART 2 OF HW WHEN X=Z

#Create matrices/variables to score key data from each simulation:
ols_gamma_hat_mtx <- matrix(0, nrow=length(true_gamma), ncol=t) #store initial OLS coefficients (used to seed optim fn)
gamma_hat_mtx <- matrix(0, nrow=length(true_gamma), ncol=t)
biased_beta_hat_mtx <- matrix(0, nrow=length(x_obs[1,]), ncol=t) #will hold naive biased OLS estimates of uncensored data
beta_hat_mtx <- matrix(0, nrow=length(x_obs[1,]), ncol=t) #will hold 2-step corrected OLS estimates of uncensored data
mu_hat_mtx <- matrix(0, nrow=1, ncol=t) #coefficient on Inverse Mills Ratio in 2nd step regression
avg_inv_mills <- 0   #used to verify calculation of Inverse Mills Ratio (initialize to zero, running sum of avg from each trial)
error_cov <- matrix(0, nrow=1, ncol=t) #store cov(u,v) each trial to take average and compare to mu_hat
avg_censored <- 0 #track number of censored data points throughout the simulation.

############## MAIN LOOP FOR 2-STEP SIMULATION (FIRST PROBIT, THEN OLS) ####################

for (i in 1:t){
  
  ########################## STEP 1: GENERATE DATA FOR SIMULATION ##########################
  
  #Generate jointly normally distributed errors for this simulation (u, v), both ~N(0,1)
  v <- matrix(data=rnorm(n, mean=0, sd=1), nrow=n, ncol=1)
  u <- 0.5*v + matrix(data=rnorm(n, mean=0, sd=1), nrow=n, ncol=1) #u = 0.5v + (zero mean noise)
  error_cov[1,i] <- cov(u,v) #store cov(u,v) to average later and compare with avg of mu_hat in stage two
  
  #Generate d's and use 'Indicator Function' via ceiling command:
  d_star <- z_data %*% true_gamma + v #latent process to determine d*
  d_indicator <- ifelse(d_star > 0, 1, 0) #assigns 1 if d*>0, 0 if d*<=0 (indicator fn)
  
  avg_censored <- avg_censored + n-sum(d_indicator) #keep running sum of censored data
  
  #Generate y_star (using true model) and y_obs (using indicator function)
  y_star <- x_obs %*% true_beta + u #latent process to determine y_star
  y_obs <- matrix(999, nrow=n, ncol=1)
  for (j in 1:n){
    y_obs[j] <- d_indicator[j] * y_star[j] #indicator determines whether y_star is observed 
  }
  
  ###########################################################################################
  
  ###################### STEP 2: RUN PROBIT REGRESSION TO GET gamma_hat #####################
  
  #First generate OLS estimate of gamma_hat to use as optim() seed: inv(z'z)z'd
  ols_gamma_hat_mtx[,i] <- solve(t(z_data) %*% z_data) %*% t(z_data) %*% d_indicator 
  
  #Do Probit to get estimate of gamma_hat
  optimization <- optim(par = ols_gamma_hat_mtx[,i], fn = Probit_LL, y = d_indicator, x = z_data, 
                        method = "Nelder-Mead", hessian = FALSE, control = list(reltol = 1e-9))
  
  gamma_hat_mtx[,i] <- optimization$par
  
  ###########################################################################################
  
  ############# STEP 3: REDUCE DATASET TO ONLY INCLUDE UNCENSORED OBSERVATIONS ##############
  #This uncensored data will be used for the inverse mills ratio calcs and stage two OLS:
  
  #Append 'observed' indicator function to x_obs and keep only rows where indicator is 1
  x_uncensored <- cbind(d_indicator, x_obs)
  x_uncensored <- x_uncensored[x_uncensored[,1]==1,] #Keeps all rows where first column equals 1
  x_uncensored <- x_uncensored[,-1] #delete first column (indicator fn no longer needed)
  
  #Repeat process with z_data to removed censored rows (needed for inverse mills ratio calcs in next step)
  z_uncensored <- cbind(d_indicator, z_data)
  z_uncensored <- z_uncensored[z_uncensored[,1]==1,] #Keeps all rows where first column equals 1
  z_uncensored <- z_uncensored[,-1] #delete first column (indicator fn no longer needed)
  
  #Similarly, keep all the rows from y_obs where first column is NOT equal to zero (i.e., keep all observed rows)
  y_uncensored <- matrix(y_obs[y_obs[,1]!=0,], nrow=sum(d_indicator), ncol=1)
  
  ###########################################################################################
  
  #### STEP 4: CALC/STORE INVERSE MILLS RATIO FROM PROBIT REGRESSION FOR UNCENSORED DATA ####
  
  #Calculate inv mills ratio for all uncensored observations at once:
  inv_mills <- phi(z_uncensored %*% gamma_hat_mtx[,i]) / PHI(z_uncensored %*% gamma_hat_mtx[,i])
  avg_inv_mills <- avg_inv_mills + sum(inv_mills)/n #keep running sum of averages of inv_mills for each trial
  
  #cat("Average Inverse Mills Ratio for trial", i, "is", sum(inv_mills)/n, "\n")

  ###########################################################################################
  
  
  ########### STEP 5: RUN OLS ON UNCENSORED DATA TO ESTIMATE beta_hat AND mu_hat ############
  
  #Run naive (biased) OLS on uncensored observations WITHOUT using Inverse Mills Ratio correction:
  biased_beta_hat_mtx[,i] <- solve(t(x_uncensored) %*% x_uncensored) %*% t(x_uncensored) %*% y_uncensored
  
  #Create new data set by appending Inverse Mills Ratio data to x_obs
  stage_two_data <- cbind(x_uncensored, inv_mills)
  
  #Run OLS to get final estimates and store values in 
  stage_two_ols <- solve(t(stage_two_data) %*% stage_two_data) %*% t(stage_two_data) %*% y_uncensored
  beta_hat_mtx[,i] <- stage_two_ols[1:length(x_obs[1,]),]
  mu_hat_mtx[,i] <- stage_two_ols[length(x_obs[1,])+1,]
  
  ###########################################################################################
}

#Output statements to summarize simulation:
cat(" SIMULATION RESULTS (n =",n,"total observations per trial, some of which were censored):\n",
    "Average gamma_hat probit values after",t,"trials:",rowMeans(gamma_hat_mtx),"\n",
    "Average (biased) beta_hat OLS values after",t,"trials:",rowMeans(biased_beta_hat_mtx),"\n",
    "Average (unbiased) beta_hat OLS values after",t,"trials:",rowMeans(beta_hat_mtx),"(using Heckman correction)\n",
    "Average mu_hat OLS value (coefficient of Inverse Mills Ratio) after",t,"trials:",rowMeans(mu_hat_mtx),"\n",
    "Average cov(u,v) (i.e., theoretical mu) value after",t,"trials:",rowMeans(error_cov),"\n",
    "Ratio of mu_hat to average theoretical mu (ratio should be close to 1):", rowMeans(mu_hat_mtx)/rowMeans(error_cov),"\n",
    "Ratio of unbiased beta_hat to true_beta (ratio should be close to 1):", rowMeans(beta_hat_mtx)/true_beta,"\n",
    "Average Inverse Mills Ratio (for uncensored observations) after", t, "trials:", avg_inv_mills/t, "\n",
    "Average number of censored data points per trial:",avg_censored/t,
    "(or",avg_censored/t/n*100,"percent of observations per trial)\n")