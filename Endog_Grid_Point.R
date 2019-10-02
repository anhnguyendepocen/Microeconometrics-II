# Shai Knight-Winnig 2017

rm(list = ls()) #clear environment 

#Import functions for Probit operations:
source('Probit_Fns.R') 

#Initialize variables and generate data to be used in all trials:
n <- 200 #Number of observations in 'sample'
t <- 500 #Number of trials in simulation

true_theta <- c(1.25)
true_beta <- c(.3, .8)
true_alpha <- c(.75)
true_delta <- c(.6, .4)

x_obs <- cbind(matrix(1, nrow=n, ncol=1), matrix(data=rnorm(n*(length(true_beta)-1), mean=0, sd=1), nrow=n, ncol=length(true_beta)-1))
z_data <- matrix(data=rnorm(n*length(true_theta), mean=0, sd=1), nrow=n, ncol=length(true_theta))


error_cov <- matrix(0, nrow=1, ncol=t) #store cov(u,v) each trial to take average and compare to mu_hat
      
for (i in 1:t){
  
  ########################## STEP 1: GENERATE DATA FOR SIMULATION ##########################
  
  #Generate jointly normally distributed errors for this simulation (u, v), both ~N(0,1)
  v <- matrix(data=rnorm(n, mean=0, sd=1), nrow=n, ncol=1)
  u <- 0.5*v + matrix(data=rnorm(n, mean=0, sd=1), nrow=n, ncol=1) #u = 0.5v + (zero mean noise)
  error_cov[1,i] <- cov(u,v) #store cov(u,v) to average later and compare with avg of mu_hat in stage two
  
  y2 <- x_obs %*% true_delta + true_theta * z_data + u          
  y1 <- x_obs %*% true_beta + true_alpha * y2 + v      
  
  y_star <- matrix(999, nrow=n, ncol=1)
  for (j in 1:n){
    y_star <- ifelse(y2 > y2[j], 1, 0)
    
    
    #First generate OLS estimate of gamma_hat to use as optim() seed: inv(z'z)z'd
    #ols_guess_matrix[,i] <- solve(t(z_data) %*% z_data) %*% t(z_data) %*% d_indicator 
    
    #Do Probit to get estimate of gamma_hat
    #optimization <- optim(par = ols_gamma_hat_mtx[,i], fn = Probit_LL, y = d_indicator, x = z_data, 
    #                      method = "Nelder-Mead", hessian = FALSE, control = list(reltol = 1e-9))
    
    rm(y_star)
  }
                 
}                 
                 
                 
                 
                 
                 
                 