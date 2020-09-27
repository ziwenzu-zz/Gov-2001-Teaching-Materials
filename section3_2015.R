##################################################################
##################################################################
# Gov 2001/1002/E-2001 Section 3
# 2-11-15
# Soledad Artiz Prillaman
##################################################################
##################################################################

# Clear our workspace
rm(list=ls())


####################################
####### Likelihood Inference #######
####################################

####################################
## PLOTTING THE LIKELIHOOD

# dexp provides the PDF of the exponential distribution
# If we set x to be 12 and lambda to be .25,
# then p(y | lambda) =
dexp(x = 12, rate = .25, log=FALSE)

# To plot the likelihood, we can plot the PDF for different
#  values of lambda
# The curve function plots the inputted function 
#	with respect to what you label as x.
# 	This means that the function you put in will be on the y-axis
#	and what you specify as x will be on the x-axis.
curve(dexp(12, rate = x), 
      xlab = expression(lambda),
      ylab = expression(paste("L(", lambda, "|y = 12)", sep = "")),
      bty = "l",
      main = "Likelihood density function",
      xlim = c(0,1))


# OR we can write a function for L( lambda | y) 
#  The inputs are lambda and y (our data)
likelihood <- function(lambda, y){
  lambda * exp(-lambda * y)
}

# We can use curve to plot this likelihood function
curve(likelihood(lambda = x, y = 12),
      xlab = expression(lambda),
      ylab = expression(paste("L(", lambda, "|y = 12)", sep = "")),
      bty = "l",
      main = "Likelihood density function",
      xlim = c(0,1))


####################################
######## Bayesian Inference ########
####################################

####################################
## PLOTTING THE PRIOR

# dgamma provides the PDF of the gamma distribution
# To plot the prior, we can plot the PDF for different
#  values of lambda
## Plot our Prior with alpha = 1 and beta = 2
curve(dgamma(x, shape = 1, rate = 2), 
      xlab = expression(lambda),
      ylab = "Density",
      bty = "l",
      main = expression(paste(alpha, "= 1 and ", beta, "= 2", sep = "")),
      xlim = c(0,5))

## Plot our Prior with alpha = 3 and beta = 2
curve(dgamma(x, shape = 3, rate = 2), 
      xlab = expression(lambda),
      ylab = "Density",
      bty = "l",
      main = expression(paste(alpha, "= 3 and ", beta, "= 2", sep = "")),
      xlim = c(0,5))

####################################
## PLOTTING THE POSTERIOR

# dgamma provides the PDF of the gamma distribution
# To plot the posterior, we can plot the PDF for different
#  values of lambda
# Plot our Posterior with alpha = 1 and beta = 2
#   ==> Gamma(2, 14)
curve(dgamma(x, shape = 2, rate = 14), 
      xlab = expression(lambda),
      ylab = "Density",
      bty = "l",
      main = "Posterior Density",
      xlim = c(0,1))


# OR we can write a function for p( lambda | y) 
#  The inputs are lambda, y (our data), beta, and alpha 
posterior.propto <- function(lambda, y, beta, alpha){
  lambda^alpha * exp(-lambda * (beta + y))
}

# We can use curve to plot this likelihood function
# Let's set alpha = 1 and beta = 2
# Note that this is not our PDF because we have ommitted 
#   the constants
curve(posterior.propto(lambda = x, y = 12, alpha = 1, beta = 2),
      xlim = c(0,1),
      xlab = expression(lambda),
      ylab = "Density",
      main = "Posterior Density",
      bty = "l")


####################################
## SUMMARIZING THE POSTERIOR

## Mean of Posterior

# 1. We can simulate to get the mean of the posterior
# Again let's assume that alpha is 1 and beta is 2
# Draw 10000 values from a Gamma(2, 14) distribution
# then calculate the mean
sim.posterior <- rgamma(10000, shape = 2, rate = 14)
mean(sim.posterior)

# Or we can calculate this analytically because the mean of a gamma is alpha / beta
mean <- 2/14


## Probability that lambda is between .2 and .4

# 1. We can integrate our proportional function and divide this by the total area from our proportional function
# We have to divide by the total density of our function because of the ommitted constants
density1 <- integrate(posterior.propto,
                     alpha = 1,
                     beta = 2,
                     y = 12,
                     lower = .2,
                     upper = .4)

total.density <- integrate(posterior.propto,
                           alpha = 1,
                           beta = 2,
                           y = 12,
                           lower = 0,
                           upper = 100)

density1$value / total.density$value # 0.2067


# 2. Or we can calculate this directly because we know the functional form of the posterior:
pgamma(.4, shape = 2, rate = 14) - pgamma(.2, shape = 2, rate = 14) # 0.2067



####################################
########## Neyman-Pearson ##########
####################################

####################################
## CALCULATING P-VALUES

# create a population with mean 12, sd 4
pop <- rnorm(n=1e6, mean=12, sd=4)

# draw a single sample of size 100 from population
my.sample <- sample(pop, size=100, replace=T)

# calculate our test statistic
# this is a test statistic for the sample mean 
#  as an estimator of the true mean (which we set to be 12)
test.statistic <- mean(my.sample) /(sd(my.sample)/10)

# find the p-value
p.value <- 2*(1-pnorm(test.statistic))
p.value