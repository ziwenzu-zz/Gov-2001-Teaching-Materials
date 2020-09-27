##################################################################
##################################################################
# Gov 2001/1002/E-2001 Section 1
# 1-28-15
# Soledad Artiz Prillaman
##################################################################
##################################################################

# Clear our workspace
rm(list=ls())


####################################
############ Simulation ############
####################################

######
# Simulating from distributions

## R has a bunch of canned functions which make this easy

## Say we want to look at the distribution of salaries
## and we assume that salaries are normally distributed 
## with mean 40000 and standard deviation 10000

## It may be hard to conceptualize this distribution,
## So we can sample from it using the rnorm() command
##    n is the sample size
salaries <- rnorm(n=1000000 , mean=40000 , sd=10000)

## Plot the distribution
plot(density(salaries))

## Now, what if we wanted to know the density of the distribution
## i.e. the height of the distribution
## when x = 20000
## dnorm provides the PDF of the normal distribution
dnorm(20000, mean=40000, sd=10000)

## Or maybe we want to know the probability that a salary is 
## is less than 20000
## pnorm provides the CDF of the normal distribution
pnorm(20000, mean=40000, sd=10000)

## Or maybe we want to know what salary marks the 95% percentile
## qnorm provides the inverse CDF of the normal distribution
qnorm(0.95, mean=40000, sd=10000)
## this tells us that under our assumed model, 
## 95% of people earn a salary of less than roughly $56000


## There are many other distributions which you can sample from in a similar way!



####################################
## EXPECTATION VIA SIMULATION

# Ex- You are waiting for the red line. How long with it take
#     for the T to get to you?
# X can be {1,2,3,4,...} minutes
# X is a geometric random variable measuring the number of minutes 
#   until the train gets to you (including the minute it gets to you)
# Assume that the probability that the T comes in any given minute is
#   pi = .2


# Now say we want to know the expected value of X - i.e. the expected 
#   amount of time you are going to wait for the train
# We can use simulation to approximate this value

# Set the seed so that you can replicate your results
#   doesn't matter what it is set to
set.seed(02138)
# draw 100000 observations from a geometric distribution
draws <- rgeom(n = 100000, prob = .2)
# and then find the mean of these to approximate the expected value
mean(draws)


# NOTE: For some versions of R, rgeom does not include the minute you arrive
#   i.e. rather than being the number of failures (minutes the train doesn't come)
#        and the first success it is the number of failures before the first success
# So we have to add 1 to rgeom to get the correct distribution - each value it gives us 
#       is the number of minutes the train did not show up and then we add 1 for the minute
#       it did show up
set.seed(02138)
# draw 100000 observations from a geometric distribution
draws <- 1 + rgeom(n = 100000, prob = .2)
# and then find the mean of these to approximate the expected value
mean(draws)


# What if instead you wanted to know the expected value of a function
# Say E[sqrt(1 + X)]
# We can calculate this by simulating X, running it through the function g(x) and finding the mean 
# We can do this because of LOTUS which says that E[g(X)] = \sum_x g(x)P(X=x)

# 1. Draw X from a geometric
set.seed(02138)
draws <- 1 + rgeom(n = 100000, prob = .2)
# 2. Calculate g(x) for each draw
g.draws <- sqrt(1 + draws)
# Find the mean of this
mean(g.draws)



#######
## Ex - Lets assume we have the same question but we want to model time as continuous,
#       i.e. the train can arrive at any point, not just on exact minutes
# X can be any value greater than 0
# X is an exponential random variable
# Assume that the rate at which X occurs is .25, i.e. lambda = .25

# Again we want to find the expected value of X using simulation
# Set the seed
set.seed(2001)
# draw 100000 observations from an exponential distribution
draws <- rexp(n = 100000, rate = .2)
# and then find the mean of these to approximate the expected value
mean(draws)


# Again we want E[sqrt(1 + X)]
# 1. Draw X from an exponential
set.seed(2001)
draws <- rexp(n = 100000, rate = .2)
# 2. Calculate g(x) for each draw
g.draws <- sqrt(1 + draws)
# Find the mean of this
mean(g.draws)



####################################
## PROBABILITY VIA SIMULATION

# Ex- We have an urn with five red balls and five green balls
# What is the probability of drawing four balls of the same color?

# 1. Construct our population
urn <- c("G","G","G","G","G","R","R","R","R","R")
urn
# Or
urn <- c(rep("red",5),rep("green",5))
urn
# Or
urn <- c(rep(1,5),rep(0,5))
urn

# 2. Take one sample
set.seed(1217)
# sample 4 elements from the vector urn without replacement
draw <- sample(x = urn, size = 4, replace = FALSE)
draw

# 3. Determine if a success or failure
sum(draw) == 4 
sum(draw) == 0
# specify that it is a success if the vector adds to 1 or 0
success <- (sum(draw) == 4) | (sum(draw) == 0)
success

# 4. Repeat with a for-loop
set.seed(1217)
# Specify the number of times you want to repeat
sims <- 100000
# Create an empty vector to store your output
success <- NULL
# Repeat steps 2-3 over and over and store your results from each repetition
#   in the vector success
for(i in 1:sims){
  draw <- sample(x = urn, size = 4, replace = FALSE)
  success[i] <- sum(draw) == 4 | sum(draw) == 0
}
head(success)

# (4. Repeat with replicate)
# Write a function that takes in a vector of draws (our sample of 4 balls)
#   and returns whether it was a success or not
eval <- function(draw){
  success <- sum(draw) == 4 | sum(draw) == 0
  return(success)
}
set.seed(1217)
sims <- 10000000
# Replicate teh function over and over and store in a new vector 
success <- replicate(sims, eval(sample(urn, 4, replace=F)))
head(success)
# Exactly the same as our for-loop! (Note: this is because we set it to the same seed)

# 5. Determine proportion of success
sum(success)/sims
# Or
mean(success)




####################################
############ FUNCTIONS #############
####################################


######### QUIZ ###########
## Write a function that takes in a vector x and a vector y 
#  and reports the r^2 from a regression of x on y
#  do not use the lm() command at all! 



















# This function takes as its inputs a vector of observations for the single covariate x
# and a vector of observations for the y
# It outputs a single value of the r2
r2.function <- function(x,y){
  # Find our estimate of b1 using the first normal equation
  b1 <- sum((x-mean(x))*y)/sum((x-mean(x))^2)
  # Find our estimate of b0 using the second normal equation
  b0 <- mean(y) - b1*mean(x)
  # using our estimates of b1 and b0 find our yhat
  yhat <- b0 + b1*x
  # use these in the formula for r2
  r2 <- sum((yhat-mean(y))^2)/sum((y-mean(y))^2)
  return(r2)
}
