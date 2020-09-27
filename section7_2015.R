##################################################################
##################################################################
# Gov 2001/1002/Stat E-200 Section 7
# 3-1-15
# Soledad Artiz Prillaman
##################################################################
##################################################################

# Clear our workspace
rm(list=ls())

# Set Working Directory
setwd("/Users/Sole/Google Drive/Gov 2001/Spring2015/Section/Section 7 - Probit and QoI (3-11)/Code/")


######################################
########## The Probit Model ##########
######################################

# Load in the Data
require(foreign)
votes <- read.dta("votes0408.dta")

## Look at the dataset
head(votes)
apply(votes, 2, summary)

# incwin: dependent variable. whether incumbent/incumbent's party won the election (1) or not (0)
# open: is the seat open (1) or is there an incumbent running (0)
# freshman: is the incumbent a freshman (1) or not (0)
# incpres: percentage of votes received by incumbent party in the district in last presidential election

# Specify our data for the model
y <- votes$incwin
X <- as.matrix(votes[,2:4])

######################################

# Write a function for the likelihood of the probit
# L(beta) = sum_i^n y*log(Phi(XB)) + (1-y)*Log(1-Phi(XB))

ll.probit <- function(beta, y=y, X=X){
  # Add a column of 1's if X doesn't already include it
  if(sum(X[,1]) != nrow(X)) X <- cbind(1,X)
  # Specify the systematic component using pnorm
  # Pnorm gives the cdf of the normal distribution, which is Phi
  # Set log=T to take the log of the CDF
  phi <- pnorm(X%*%beta, log = TRUE) 
  # Calculate 1-Phi by setting lower.tail = F
  opp.phi <- pnorm(X%*%beta, log = TRUE, lower.tail = FALSE)
  # Calculate the log likelihood
  logl <- sum(y*phi + (1-y)*opp.phi)
  return(logl)
}



## Test that the function works:
ll.probit(beta = rep(0, ncol(X) + 1),
         X = X,
         y = y) 
# -730.5771
# This is the same as we found with the logit!


## Estimate the MLE:
opt <- optim(par = rep(0, ncol(X) + 1),
             fn = ll.probit,
             y = y,
             X = X,
             control = list(fnscale = -1),
             hessian = T,
             method = "BFGS")

# Extract the coefficients from optim
coefs <- opt$par 
# Calculate the variance
varcov <- -solve(opt$hessian)
# Calculate the standard errors
ses <- sqrt(diag(varcov))

# Format into a nice table
res <- cbind(coefs, ses)
colnames(res) <- c("Coefficient", "SE")
rownames(res) <- c("Intercept", "Open Seat", "Freshman", "Pres Vote")

require(xtable)
xtable(res, caption="DV: Incumbent Election Win", digits=4)

######################################
####### Quantities of Interest #######
######################################

#1. We have our estimates of Betahat and its variance

#2. Draw our beta tildes from a multivariate normal distribution
#     with mean and variance = our estimates
#     to account for estimation uncertainty
#   Note: make sure to use the entire variance covariance matrix!
library(mvtnorm)
set.seed(1234)
sim.betas <- rmvnorm(n = 10000,
                     mean = coefs,
                     sigma = varcov) 
head(sim.betas)
dim(sim.betas) ## 10,000 by 4 matrix of simulated betas

#3. Set our values of X = Xc
# Let's identify a case we care about
# Say when people are up for reelection
# They were a freshmen
# And their presidential party got 48% of the vote
# Don't forget about the intercept!
Xc <- c(1,1,1,48)

# Alternatively, we could just set every value at its median
Xc <- apply(X = cbind(1,X), MARGIN = 2, FUN = median)


######################################
# PREDICTED VALUES

#4. Calculate the pi. systematic component
# and
#5. Draw our Y tildes from the stochastic component
#     to account for fundamental uncertainty
pv.ests <- c() # Create an empty vector
for(i in 1:10000){ 
  pi.tilde <- pnorm(Xc%*%sim.betas[i,])
  pv.ests[i] <- rbinom(1, 1, pi.tilde)
}
head(pv.ests)

# Look at the results
hist(pv.ests,
     main = "Predicted Values",
     xlab = "Predicted Value of Y")


######################################
# EXPECTED VALUES

#4. Calculate the pi. systematic component
# and
#5. Draw our Y tildes from the stochastic component
#     to account for fundamental uncertainty
#6. Calculate the expected value of Y 
ev.ests <- c() # Create an empty vector
for(i in 1:10000){ 
  pi.tilde <- pnorm(Xc%*%sim.betas[i,])
  y.ests <- rbinom(10000, 1, pi.tilde)
  ev.ests[i] <- mean(y.ests)
}
head(ev.ests)

# Look at the results
hist(ev.ests,
     main = "Expected Values",
     xlab = "Predicted Probability of Incumbent Win")

mean(ev.ests)
quantile(ev.ests, c(.025,.975))


# We can use a shortcut though
set.seed(1234)
sim.betas <- rmvnorm(n = 10000,
                     mean = coefs,
                     sigma = varcov) 
ev.ests2 <- pnorm(Xc%*%t(sim.betas))  
# Because expected values average over the fundamental uncertainty
# we don't have to include the simulation from the stochastic component
# IF E[Y] = theta (the parameter)
mean(ev.ests2)
quantile(ev.ests2, c(.025,.975))


#What if I wanted to calculate the expected values for
# many values of the covariates?
#Let's do this for president vote shares from 30% to 70%
presvals <- 30:70
# Create an empty matrix so that we can put out expected values for each level of presvals
# Remember we get a 10000 element vector for each
ev.ests <- matrix(data = NA, ncol = length(presvals),
                 nrow=10000)

# Loop over all values of presvals
for(j in 1:length(presvals)){
  Xc.new <- Xc
  Xc.new[4] <- presvals[j]
  ev.ests[,j] <- pnorm(Xc.new%*%t(sim.betas)) # Notice we don't have to simulate new betas
}

plot(presvals, apply(ev.ests,2,mean), ylim=c(0,1), 
     xlab="Presidential Vote Share", ylab="Probability of Incumbent Win", main="Expected Values")
segments(x0 = presvals, x1 = presvals, 
         y0 = apply(ev.ests, 2, quantile, .025), 
         y1 = apply(ev.ests, 2, quantile, .975))



######################################
# FIRST DIFFERENCES

# Let's compare the expected value for freshmen and non
# Xc already is set for freshmen
Xc.fresh <- Xc
Xc.nofresh <- Xc
Xc.nofresh[3] <- 0

# Calculate the first differences as the 
# difference in the expected values
# Remember we can only do this because the expected value
# is pi
fd.ests <- pnorm(Xc.fresh%*%t(sim.betas)) -
  pnorm(Xc.nofresh%*%t(sim.betas))

# Look at the results
hist(fd.ests,
     main = "First Differences",
     xlab = "Difference in Predicted Probability \n for Freshmen and Non-freshmen")

mean(fd.ests)
quantile(fd.ests, c(.025,.975))

######################################
# USING ZELIG

# We are going to be using Zelig 5.0, which is still in beta release
# Direct any generic questions you have about the software to:
#   https://groups.google.com/forum/#!forum/zelig-statistical-software
# Report any bugs you find to:
#   https://github.com/IQSS/Zelig/issues

# Make sure you have all of the necessary dependencies installed
install.packages("jsonlite")
install.packages("AER")
install.packages("dplyr")
install.packages("quantreg")
install.packages("geepack")
install.packages("maxLik")
# To install Zelig 5.0, you must type exactly:
install.packages("Zelig", type = "source", repos = "http://r.iq.harvard.edu/")
library(Zelig)
# If you want more detailed instructions, go to http://docs.zeligproject.org/en/latest/installation_quickstart.html


# Now let's run our model
# This is very similar to the lm function
# but you specify the model as probit
zelig.out <- zelig(incwin ~ open + freshman + incpres, 
                   data = votes, model = "probit")
summary(zelig.out)

# Set the values of the covariates
x.evs <- setx(zelig.out, open = 1, freshman = 1, incpres = 48)
x.evs

# Simulate predicted probabilities
set.seed(1234)
zelig.sim <- sim(zelig.out, x=x.evs)
summary(zelig.sim)
# Plot our simulations
par(mar=c(2.5,2.5,2.5,2.5)) # Set our plot margins
plot(zelig.sim)


######################################
########## Common Questions ##########
######################################

######################################
# Comparing probit and logit

# Run our model with a logit rather than a probit systematic component
ll.logit <- function(beta, y, X){
  if(!all(X[,1] == 1)){X <- cbind(1,X)}
  xb <- X %*% beta
  logl <- sum(y * xb - xb - log(1 + exp(-xb)))
  return(logl)
}
opt.l <- optim(par = rep(0, ncol(X) + 1),
               fn = ll.logit,
               y = y,
               X = X,
               control = list(fnscale = -1),
               hessian = T,
               method = "BFGS")
coefs.l <- opt.l$par
varcov.l <- -solve(opt.l$hessian)
ses.l <- sqrt(diag(varcov.l))


# Let's compare the coefficients
lp.coefs <- cbind(coefs, coefs.l)
colnames(lp.coefs) <- c("Probit", "Logit")
rownames(lp.coefs) <- c("Intercept", "Open Seat", "Freshman", "Pres Vote")
xtable(lp.coefs, digits=4)

# Let's compare quantities of interest (see below)
# Difference 45% presidential vote share and 55%
Xc.lo <- c(1,1,1,45)
Xc.hi <- c(1,1,1,55) 
fd.ests.p <- pnorm(Xc.hi%*%t(sim.betas)) - pnorm(Xc.lo%*%t(sim.betas))
set.seed(1234)
sim.betas.l <- rmvnorm(n = 10000, mean = coefs.l, sigma = varcov.l) 
fd.ests.l <- pnorm(Xc.hi%*%t(sim.betas.l)) - pnorm(Xc.lo%*%t(sim.betas.l))
lp.fds <- cbind(c(mean(fd.ests.p),quantile(fd.ests.p, c(.025,.975))) , c(mean(fd.ests.l),quantile(fd.ests.l, c(.025,.975))))
colnames(lp.fds) <- c("Probit", "Logit")
rownames(lp.fds) <- c("Mean", "2.5% Quantile", "97.5% Quantile")
xtable(lp.fds, digits=4, caption = "First Differences")

# What about for the tails?
# Difference 30% presidential vote share and 90%
Xc.lo <- c(1,1,1,30)
Xc.hi <- c(1,1,1,90) 
fd.ests.p <- pnorm(Xc.hi%*%t(sim.betas)) - pnorm(Xc.lo%*%t(sim.betas))
set.seed(1234)
sim.betas.l <- rmvnorm(n = 10000, mean = coefs.l, sigma = varcov.l) 
fd.ests.l <- pnorm(Xc.hi%*%t(sim.betas.l)) - pnorm(Xc.lo%*%t(sim.betas.l))
lp.fds <- cbind(c(mean(fd.ests.p),quantile(fd.ests.p, c(.025,.975))) , c(mean(fd.ests.l),quantile(fd.ests.l, c(.025,.975))))
colnames(lp.fds) <- c("Probit", "Logit")
rownames(lp.fds) <- c("Mean", "2.5% Quantile", "97.5% Quantile")
xtable(lp.fds, digits=4, caption = "First Differences")


######################################
####### Robust Standard Errors #######
######################################

# Load some data from the Gov 2001 Code library
load("Gov2001CodeLibrary.RData")
#This particular dataset is from the paper 
# "When preferences and commitments collide: 
# The effect of relative partisan shifts on International treaty compliance." 
# in International Organization by Joseph Grieco, Christopher Gelpi, and Camber Warren.

# Replicate their model
fmla <- as.formula(restrict ~ art8 + shift_left + flexible +
                     gnpcap + regnorm + gdpgrow + resgdp + bopgdp + 
                     useimfcr + surveil + univers + resvol + totvol + 
                     tradedep + military + termlim + parli + lastrest + 
                     lastrest2 + lastrest3)

# The GLM is a function which performs our maximum likelihood
# estimation
fit <-glm(fmla, data=treaty1, 
          family=binomial(link="logit")) 

# The bread is just the variance covariance matrix
# We can extract this with the vcov function 
bread <-vcov(fit)

# We can calculate the score using the estfun function
# in the sandwich package
library(sandwich)
est.fun <- estfun(fit)

#The mean is the score twice
meat <- t(est.fun)%*%est.fun
# the sandwich, which has our robust standard errors
# is the bread * meat * bread
sandwich <- bread%*%meat%*%bread
dim(sandwich)

# And now we can put our robust standard errors in our table
library(lmtest)
coeftest(fit, sandwich)

# For linear models, est fun behaves a little off
# You can use
robust <- sandwich(lm.1, meat=crossprod(est.fun)/N)


######################################
# Cluster-Robust

# First, specify the clusters
fc <- treaty1$imf_ccode
m <- length(unique(fc))
k <- length(coef(fit))

# now we have to sum up the meat by the cluster
u <- estfun(fit)
u.clust <- matrix(NA, nrow=m, ncol=k)
for(j in 1:k){
  u.clust[,j] <- tapply(u[,j], fc, sum)
}
dim(u) # n x k
dim(u.clust) # m x k

# Calculate the cluster-robust variance covariance matrix
# bread * cluster.meat * bread
cl.vcov <- bread %*% ((m / (m-1)) * t(u.clust) %*% (u.clust)) %*% bread

# And add to our table
coeftest(fit, cl.vcov)







