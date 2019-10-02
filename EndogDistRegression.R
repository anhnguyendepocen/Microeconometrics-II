# Shai Knight-Winnig 2017
# Endogenous Conditional CDF Estimation via  Monte Carlo simulation

rm(list = ls()) #clear environment 

#Start timing of code:
start_time <- Sys.time()

#Import functions for Probit operations:
source('Probit_Fns.R') 

#Initialize variables and generate data to be used in all trials:
n <- 200 #Number of observations in 'sample'
t <- 10 #Number of trials in simulation
set.seed(21)  #Uncomment this line if you want to use the same seeded random sequence as someone else to compare results

true_theta <- c(1.25) 
true_beta <- c(0.3, 0.8) 
true_alpha <- c(0.75)
true_delta <- c(0.6, 0.4)

#Create x data to be used in each simulation's 2nd step regression (2 regressors plus an intercept)
x_obs <- cbind(matrix(1, nrow=n, ncol=1), matrix(data=rnorm(n*(length(true_beta)-1), mean=0, sd=1), 
                                                 nrow=n, ncol=length(true_beta)-1))

#Create z data to be used in each simulation's 1st step Probit (x_obs plus one more regressor)
z_data <- matrix(data=rnorm(n*length(true_theta), mean=0, sd=1), nrow=n, ncol=length(true_theta))

#Create matrices/variables to score key data from each simulation:
cond_cdf_mtx <- matrix(0, nrow=n, ncol=t)
biased_beta_mtx <- matrix(0, nrow=length(true_beta), ncol=t) #will hold naive biased OLS estimates endogenous
biased_alpha_mtx <- matrix(0, nrow=length(true_alpha), ncol=t) #coefficient on y2 in 2nd step regression
beta_hat_mtx <- matrix(0, nrow=length(true_beta), ncol=t) #will hold 2-step corrected OLS estimates
alpha_hat_mtx <- matrix(0, nrow=length(true_alpha), ncol=t) #coefficient on y2 in 2nd step regression
gamma_hat_mtx <- matrix(0, nrow=1, ncol=t) #coefficient on cond CDF in 2nd step regression
error_cov <- matrix(0, nrow=1, ncol=t) #store cov(u,v) each trial

######## MAIN LOOP FOR 2-STEP SIMULATION (FIRST COND CDF ESTIMATION, THEN OLS WITH CONTROL FUNCTION) ###########

for (i in 1:t){
  
  ########################## STEP 1: GENERATE DATA FOR SIMULATION ##########################
  
  #Generate jointly normally distributed errors for this simulation (u, v), both ~N(0,1)
  v <- matrix(data=rnorm(n, mean=0, sd=1), nrow=n, ncol=1)
  u <- 0.5*v + matrix(data=rnorm(n, mean=0, sd=1), nrow=n, ncol=1) #u = 0.5v + (zero mean noise)
  error_cov[1,i] <- cov(u,v) #store cov(u,v) to average later and compare with avg of mu_hat in stage two
  
  #Generate y1 and y2 data:
  y2 <- x_obs %*% true_delta + z_data %*% true_theta + u
  y1 <- x_obs %*% true_beta + y2 %*% true_alpha + v
  
  ###################### STEP 2: RUN PROBIT REGRESSION TO GET Coefficients for Conditional CDFs #####################
  
  #for each y2, use y2 as a grid point and run probit at that grid point
  
  #matrices for storing probit coefficients for each grid point (i.e., each y2)
  probit_delta_mtx <- matrix(0, nrow=length(true_delta), ncol=n)
  probit_theta_mtx <- matrix(0, nrow=length(true_theta), ncol=n)
  
  for (j in 1:n){
    
    y_star <- ifelse(y2 > y2[j], 1, 0) #assign 1's and 0's based on above or below grid point.
    
    #Do PROBIT to get delta and theta estimates:
    #First concatenate data and generate OLS estimate use as optim() seed: 
    probit_data <- cbind(x_obs, z_data)
    ols_guess <- solve(t(probit_data) %*% probit_data) %*% t(probit_data) %*% y_star 
    
    #Do Probit to get estimate of delta, theta
    optimization <- optim(par = ols_guess, fn = Probit_LL, y = y_star, x = probit_data, 
                          method = "Nelder-Mead", hessian = FALSE, control = list(reltol = 1e-9))
    
    probit_result <- optimization$par
    
    #Store theta and delta estimates for this grid point (below two commands accommodate different lengths for delta/theta)
    probit_delta_mtx[,j] <- probit_result[1:length(true_delta),]
    probit_theta_mtx[,j] <- probit_result[(length(true_delta)+1):(length(true_delta)+length(true_theta)),]
    
    rm(y_star, probit_data, ols_guess, probit_result) #clear for next iteration
  }
  
  ############# STEP 3: COMPUTE 'RESIDUALS' USING CONDITIONAL CDFs TO BE USED AS CONTROL FUNCTION ##############
  
  #Calculate conditional CDF values at grid points using probit coefficients from step 1:
  #This will be the "control" function values that remove endogeneity
  #For each row of data, this is calculating PHI(y2-x*delta_hat-z*theta_hat)
  for (k in 1:n){
    cond_cdf_mtx[k,i] <- PHI(y2[k] - x_obs[k,]%*%probit_delta_mtx[,k] - z_data[k,]%*%probit_theta_mtx[,k])
  }
  
  ########### STEP 4: RUN OLS ON Y1 DATA WITH CONTROL FUNCTION INCLUDED TO OBTAIN CONSISTENT ESTIMATES ############
  
  #Run naive (biased) OLS WITHOUT using control function correction:
  biased_data <- cbind(x_obs, y2)
  biased_ols <- solve(t(biased_data) %*% biased_data) %*% t(biased_data) %*% y1
  biased_beta_mtx[,i] <- biased_ols[1:length(true_beta),]
  biased_alpha_mtx[,i] <- biased_ols[length(true_beta)+1,]
  
  #Create new data set by appending Inverse Mills Ratio data to x_obs
  stage_two_data <- cbind(x_obs, y2, cond_cdf_mtx[,i])
  
  #Run OLS to get final estimates and store values
  stage_two_ols <- solve(t(stage_two_data) %*% stage_two_data) %*% t(stage_two_data) %*% y1
  beta_hat_mtx[,i] <- stage_two_ols[1:length(true_beta),]
  alpha_hat_mtx[,i] <- stage_two_ols[length(true_beta)+1,]
  gamma_hat_mtx[,i] <- stage_two_ols[length(true_beta)+2,]
  

  #clear for next trial:
  rm(probit_delta_mtx, probit_theta_mtx, biased_data, biased_ols, stage_two_data, stage_two_ols, u, v, y1, y2)
}

#Stop timing of code:
end_time <- Sys.time()
total_time <- as.numeric(end_time-start_time, units="secs")

#Output statements to summarize simulation:
cat(" SIMULATION RESULTS:\n",
    "Total Processing Time:",total_time/60,"minutes\n",
    "Number of Trials:",t,"\n",
    "Total Observations Per Trial:",n,"\n",
    "Average (biased) beta_hat OLS values after",t,"trials:",rowMeans(biased_beta_mtx),"\n",
    "Average (unbiased) beta_hat OLS values after",t,"trials:",rowMeans(beta_hat_mtx),"(using control function correction)\n",
    "Average alpha_hat value (coefficient on y2) after",t,"trials:",rowMeans(alpha_hat_mtx),"\n",
    "Average gamma_hat value (coefficient of control fn residual) after",t,"trials:",rowMeans(gamma_hat_mtx),"\n",
    "Average cov(u,v) value after",t,"trials:",rowMeans(error_cov),"\n",
    "Ratio of unbiased beta_hat to true_beta (ratio should be close to 1):", rowMeans(beta_hat_mtx)/true_beta,"\n",
    "Ratio of alpha_hat to true_alpha (ratio should be close to 1):", rowMeans(alpha_hat_mtx)/true_alpha,"\n")
