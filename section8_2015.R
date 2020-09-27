##################################################################
##################################################################
# Gov 2001/1002/ Stat E-200 Section 8
# 3-25-15
# Soledad Artiz Prillaman
##################################################################
##################################################################

# Clear our workspace
rm(list=ls())

# Set Working Directory
setwd("")

######################################
###### The Ordered Probit Model ######
######################################

# We are going to estimate an ordered probit model of 
# cooperation on economic sanctions 
# To do this with zelig, you just specify the model as "oprobit"
# For order logit, specify "ologit"
# NOTE: You may have to run the following code to run an oprobit model
install.packages("VGAM")
library(Zelig)
library(devtools)
devtools::install_github("IQSS/ZeligChoice")
library(ZeligChoice)

# Load in the sanction data from the Zelig package
data(sanction)
head(sanction)
# We see that coop looks like an ordered categorical variable


# Run the model on factor(coop)
# factor() is telling R to turn coop into a matrix (but it does not do this in the data set)
z.out <- zelig(factor(coop) ~ target + cost + mil, 
               model="oprobit", data=sanction)
# Let's look at this output
z.out$result
names(z.out)

# But we don't know what this means

# Suppose we want to compare cooperation when there is or is not military action in addition to the sanction
# We set the values of our covariates
x.low <- setx(z.out, mil = 0) 
x.high <- setx(z.out, mil = 1)

# Now let's simulate!
# set x1 to x.high and x=x.low to compare these two
s.out <- sim(z.out, x = x.low, x1 = x.high)
summary(s.out)
# This gives us expected values for our pi (ie predicted probabilities) given x.low and x.high
# and this gives us the first differences between them

# We can plot this to get a lot of cool things
plot(s.out)


######################################

## We can do this WITHOUT USING ZELIG

# First we want to ake a matrix for the y's indicating what category it is in
y <- sanction$coop 
# Find all of the unique values of y and order them from smallest to largest
y0 <- sort(unique(y))
m <- length(y0)
#Create an empty matrix
Z <- matrix(NA, nrow(sanction), m)
# Fill in our matrix with logical values if the observed value is in each category
# Remember R can treat logical values as 0/1s
for (j in 1:m){Z[,j] <- y==y0[j]}
X <- cbind(sanction$target, sanction$cost, sanction$mil)


# Then we want to write a function for the log-likelihood
ll.oprobit <- function(par, Z, X){
  beta <- par[1:ncol(X)]
  tau <- par[(ncol(X)+1):length(par)]
  ystarmu <- X%*%beta
  m <- length(tau) + 1
  probs =cprobs = matrix(nrow=length(ystarmu), ncol=m)
  for (j in 1:(m-1)){cprobs[,j] <- pnorm(tau[j]- ystarmu)}
  probs[,m] <- 1-cprobs[,m-1]
  probs[,1] <- cprobs[,1]
  for (j in 2:(m-1)){probs[,j] <- cprobs[,j] - cprobs[,(j-1)]}
  sum(log(probs[Z]))
}


# Last we optimize the log-likelihood
par <- c(rep(1,3),0,1,2)
out <- optim(par, ll.oprobit, Z=Z, X=X, 
      method="BFGS", control=list(fnscale=-1), hessian=T)

out$par
sqrt(diag(solve(-out$hessian)))

# Compare to Zelig
z.out$result





######################################
### The Zero-Inflated Logit Model ####
######################################

# Load Data
fish <- read.table("http://www.ats.ucla.edu/stat/data/fish.csv", sep=",", header=T)

# Our covariates on whether or not you catch a fish given that you are fishing are
# whether or not you have a child and how many people in your group
X <- fish[c("child", "persons")]
# Our covariates on whether or not you fish is just how many people in your group
Z <- fish[c("persons")]
# Add intercepts
X <- as.matrix(cbind(1,X))
Z <- as.matrix(cbind(1,Z))
# Make Y a dummy
y <- ifelse(fish$count>0,1,0)



# Write a function for the log-likelihood
ll.zilogit <- function(par, X, Z, y){
  beta <- par[1:ncol(X)]
  gamma <- par[(ncol(X)+1):length(par)]
  phi <- 1/(1+exp(-Z%*%gamma))
  pie <- 1/(1+exp(-X%*%beta))
  sum(y*log((1-phi)*pie) + (1-y)*(log(phi + (1-phi)*(1-pie))))
}


# Then we can optimize our likelihood function
par <- rep(1,(ncol(X)+ncol(Z)))
out <- optim(par, ll.zilogit, Z=Z, X=X,y=y, method="BFGS", 
             control=list(fnscale=-1), hessian=TRUE)

out$par


# These coefficients don't mean much, so we want to simulate
# to get the predicted probability of Not Fishing
varcv.par <- solve(-out$hessian)
library(mvtnorm)
# Simulate our gammas and our betas
sim.pars <- rmvnorm(10000, out$par, varcv.par)
# Subset to only the parameters we need (gammas)
# Better to simulate all though
sim.z <- sim.pars[,(ncol(X)+1):length(par)]


# Calculate the predicted probabilities for four different group sizes
person.vec <- seq(1,4)
# Each row is a different vector of covariates (our setx)
Zcovariates <- cbind(1, person.vec)
# Calculate predicted probability (inverse logit) for each set of covariates
# Remember that for each group size, this is a vector of 10000 predicted probabilities
exp.holder <- matrix(NA, ncol=4, nrow=10000)
for(i in 1:length(person.vec)){
  exp.holder[,i] <- 1/(1+exp(-Zcovariates[i,]%*%t(sim.z)))
}

# Let's plot the distribution for each group size
plot(density(exp.holder[,4]), col="blue", xlim=c(0,1), 
     main="Probability of a Structural Zero", xlab="Probability")
lines(density(exp.holder[,3]), col="red")
lines(density(exp.holder[,2]), col="green")
lines(density(exp.holder[,1]), col="black")
legend(.7,12, legend=c("One Person", "Two People", 
                       "Three People", "Four People"), 
       col=c("black", "green", "red", "blue"), lty=1)




######################################
######## The Binomial Model ##########
######################################

# Let's create some data
x1 <- rnorm(1000,0,1)
x2 <- rnorm(1000,9,.5)
pi <- inv.logit(-5 + .4*x1 +.6*x2)
y <- rbinom(1000,10,pi)

# Write a function for our log-likelihood
ll.binom <- function(par, N, X, y){
  pi <- 1/(1 + exp(-1*X%*%par))
  out <- sum(y * log(pi) + (N - y)*log(1-pi))
  return(out)
}

# Then we can optimize the likelihood
my.optim <- optim(par = c(0,0,0), fn = ll.binom, 
                  y = y, X = cbind(1,x1,x2), N = 10, 
                  method = "BFGS", control=list(fnscale=-1), hessian=T)

my.optim$par
# This looks good given what we set the coefficients to when we defined pi


