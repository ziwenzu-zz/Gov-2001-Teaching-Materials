##################################################################
##################################################################
# Gov 2001/1002/ Stat E-200 Section 5
# 2-25-15
# Soledad Artiz Prillaman
##################################################################
##################################################################

# Clear our workspace
rm(list=ls())


#####################################
### Maximum Likelihood Estimation ###
#####################################

## UNIVARIATE EXAMPLE

## Calculate the maximum likelihood estimator 
# Model:
#       y ~ Expo(\lambda_i)
#       \lambda_i = \lambda
# Data:
#       Y: (1,5,8,2,1)


# First, let's create a function for our log-likelihood
# This function takes as its inputs 
#     par - the parameter lambda
#     y - the data
ll.expo <- function(lambda, y){
  length(y)*log(lambda) - lambda*sum(y)
}

# Create a vector of our data (5 observations and no covariates)
Y <- c(1,5,8,2,1)


# Now, let's use optim to find the MLE
# Optim takes as its inputs
#     par - the starting values for all of the parameters
#           this is where the algorithm starts
#           if there is only one parameter as in our example, this is single value
#           if there are multiple parameters, this is a vector with values for each parameter
#     fn - the function that you want to optimize
#     y - any inputs the function takes that aren't your parameter 
#         you would specify these by their name in the function 
#         (eg. if we had put data in our function instead of y this would be data = Y)
#     control = list(fnscale = -1) -- optim automatically searches for a minimum
#           by specifying to take -1, it tricks it to find the maximum
#     method - the algorithm you want to choose
#              other options: CG, L-BFGS-B, SANN, Nelder-Mead
#     hessian = TRUE -- tells optim to calculate the hessian matrix at the mle estimate
opt.expo <- optim(par = 0.01,  
                  fn = ll.expo, 
                  y = Y,
                  control = list(fnscale = -1), 
                  method = "BFGS", 
                  hessian = TRUE)
opt.expo

# If we want to pull out our MLE
mle <- opt.expo$par
mle
# This is the same as 5/17, our analytic MLE

# If we want to pull out the matrix of second derivatives evaluated at the MLE
hessian <- opt.expo$hessian
hessian
# This is the same as -289/5, our analytic hessian


#####################################
# The CLT and MLE

# Let's convince ourselves that as we take many samples of growing size, 
# the distribution of our MLE estimates approximates a normal distribution
# with mean of the true value

# We calculated analytically that the MLE is
true <- 5/17

# Set our sample size
n <- 10

# Draw 1000 samples of size n from the exponential distribution with rate = the truth
set.seed(1234)
data <- sapply(seq(1,1000),function(x)rexp(n, rate=true))

# Write a function for our log-likelihood - this is the same as the function we wrote before
llexp <- function(param, y){sum(dexp(y, rate=param, log=T))}

# Create an empty vector to store our estimates of lambda hat
out <- NULL
# For each of our 1000 samples, calculate the mle and store
for(i in 1:1000){
  out[i] <- optim(par=c(1), fn=llexp, y=data[,i],
                  method="BFGS", control=list(fnscale=-1))$par
}

#Look at a histogram of the estimates
hist(out, 
     breaks = 30, 
     xlab = expression(lambda), 
     main = expression(paste("Histogram of ", lambda, " for n = 10", sep = "")) )

## Try this for different values of n and see what happens!



#####################################

## MULTIVARIATE EXAMPLE


## Calculate the maximum likelihood estimator 
# Model:
#       y ~ Expo(\lambda_i)
#       \lambda_i = 1/exp(X_i\beta)


# Let's create some fake data
set.seed(02139)
n <- 1000
# sample whether or not it was a Friday
Friday <- sample(c(0,1), n, replace=T)
# sample the number of minutes behind schedule
minsSch <- rnorm(n, 3, .5)
# Create Y using our covariates and what we are "deciding" are the true betas
Y <- rexp(n, rate = 1/exp(1.25 - .5*Friday +.2*minsSch))
# Combind all to get data frame
data <- as.data.frame(cbind(Y, Friday, minsSch))

# We can look at our data and see that it looks exponential
hist(Y, col = "goldenrod", main = "Distribution of y")


# Let's create a function for our log-likelihood
# This function takes as its inputs 
#     param - the parameters beta0, beta1, and beta2
#     y - a vector of the data on wait time
#     x - a matrix of the covariate data
llexp <- function(param, y, x){
  rate <- 1/exp(x%*%param)
  sum(dexp(y, rate=rate, log=T))
}
# Alternatively we could write this function like this
llexp2 <- function(param, y,x){
  cov <- x%*%param
  sum(-cov - 1/exp(cov)*y)
}


#Create our matrix of X with an intercept
X <- cbind(1, Friday, minsSch)

#Specify starting values for all of our parameters
param <- c(1,1,1)

#Solve using optim
out <- optim(param, 
             fn=llexp, 
             y=Y, 
             x=X, 
             method="BFGS", 
             hessian=T, 
             control=list(fnscale=-1))

# Get the MLE
out$par

## To get the variance:
# Get the Hessian from optim
H <- out$hessian

# Calculate the observed fisher information
I <- -H

# Calculate the variance-covariance matrix
# The variances of our parameters are in the diagonal
# and the covariances are in the off-diagonal
V <- solve(I)

# Get the standard errors
ses <- sqrt(diag(V))
ses



#####################################
####### Likelihood Ratio Test #######
#####################################

# Given the model we used in the multivariate example above
# We want to test whether or not including the minutes behind schedule
# changes the mle

# We can do this using a likelihood ratio test with
# Unrestricted model: Wait = \beta_0 + \beta_1*Minutes + \beta_2*Friday
# Restricted model: Wait = \beta_0 + \beta2*Friday

# Calculate the values of the likelihood at the MLE under the 
# UNRESTRICTED model
unrestricted <- optim(param, 
                      fn=llexp, 
                      y=Y, 
                      x=X, 
                      method="BFGS", 
                      hessian=T, 
                      control=list(fnscale=-1))
# Note that what we care about is the height of the log-likelihood fuction
# not the value of the parameter that maximizes the likelihood function
unrestricted$value


# Calculate the values of the likelihood at the MLE under the 
# RESTRICTED model
restricted <- optim(c(1,1), 
                    fn=llexp,
                    y=Y, 
                    x=cbind(1, Friday), # We exclude minutes data 
                    method="BFGS", 
                    hessian=T, 
                    control=list(fnscale=-1))
restricted$value


# Calculate our test statistic given these values of the log-likelihood
# We don't have to take the log because we stored our log-likelihood values
r <- 2*(unrestricted$value - restricted$value)
r

# Calculate the p-value for this test statistic
# we have 1 degree of freedom because we have only 1 restriction (we are exluding only one variable)
1-pchisq(r,df=1)
# We reject our null hypothesis that minutes didn't matter













