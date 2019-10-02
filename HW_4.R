# Shai Knight-Winnig 2017

rm(list = ls()) #clear environment

n<-2000
x <- cbind(matrix(rnorm(n,1,2),n,1),matrix(1,n,1))
beta <- c(-.5,1)
u <- matrix(rnorm(n, 0, 1),n,1)

xbeta <- x%*%beta

y <- xbeta + u


sx <- sd(xbeta) #std dev of our random variable 
bin_width <- 1.06 #the higher the larger the bin width and the "smoother" the function this value is optimal for Gaussian
fxn <- matrix(0,n,1)
fxmm <- matrix(0,n,1)
h <- bin_width*sx*(n^(-.2)) #"bin" size, basically the width of our kernel, this is somehow the optimal value for Gaussian
for (k in 1:n){
  
  ##This calcs the "z-score" of each x with 0 at x_i
  for(h in 1:n){
    
    d = xbeta[k]-xbeta
    
    
  }
  
  ##Builds the kernel using the std normal pdf and inputing the z score of each x
  fx1 <- dnorm(d, mean = 0, sd = 1)
  num <- fx1*y
  fxn[k] <- colMeans(num)/colMeans(fx1)  
  

}

##Check accuracy of expectation
error <- matrix(0,n,1)
for(i in 1:n){
  
  error[i] <- fxn[i]-xbeta[i]
  
}

cat("Avg error", colMeans(error)/colMeans(xbeta))