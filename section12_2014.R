##################################################################
##################################################################
# Gov 2001/1002/E-2001 Section 12
# 4-23-14
# Soledad Artiz Prillaman
##################################################################
##################################################################

# Clear our workspace
rm(list=ls())

library(Zelig)
library(systemfit)
library("ZeligChoice")

###################################################
###### Seemingly Unrelated Regression Models ######
###################################################

# Load the grunfeld data from the Zelig packages
data(grunfeld)
names(grunfeld)
summary(grunfeld)


# Run the Seemingly Unrelated Regression for GE investment and 
# Westinghouse investment
# Store our two equations then run model
# To run this model we can use the systemfit package
formula <- list(mu1 = Ige ~ Fge + Cge, 
                mu2 = Iw ~ Fw + Cw)
sur.out <- systemfit(formula = formula, 
               method = "SUR", data = grunfeld)
sur.out

# To get the variance - covariance matrix 
vcov.sur<-vcov(sur.out)
vcov.sur
# We can see that there are dependencies in our parameters across models



#######################################
###### Multinomial Choice Models ######
#######################################

# Load the mexico data from the Zelig package
# This has data on vote choice across three parties in Mexico
data(mexico)
names(mexico)
summary(mexico$vote88)


# Let's run a multinomial logit model of voting behavior
# Our covariates will be age and gender
# Since vote choice is a categorical variable, we need to factorize it
# to get a matrix of y's
ml.out <- zelig(as.factor(vote88) ~ age + female, 
               model = "mlogit", data = mexico)
ml.out$result

#Let's look at the difference in expected vote choice for the youngest and the oldest
x.young <- setx(ml.out, age = min(mexico$age)) 
x.old<- setx(ml.out, age = max(mexico$age))
ml.sim <- sim(ml.out, x1 = x.old, x = x.young)
summary(ml.sim)

# Let's look at our results
plot(ml.sim)



##########################
###### Missing Data ######
##########################

# Let's load the Amelia package and the africa data in this package
library(Amelia)
data(africa)
names(africa)
summary(africa) # we have a couple of missing values in gdp pc and trade

africa[1,]
year      country gdp_pc  infl trade civlib population
1 1972 Burkina Faso    377 -2.92 29.69    0.5    5848380

# Let's run our imputation
# x specifies the data to be imputed
# all variables in x will be used as predictors in the imputation model
# cs and ts specify if cross-sectional or time-series
# logs specifies variables which need to be logged
# m specifies the number of imputed data sets
# the greater m is, the less likely you are to have outliers
set.seed(1234)
a.out <- amelia(x = africa, cs = "country", ts = "year", 
                logs = "gdp_pc", m = 5)

# If you wanted to save these 5 data sets, you could:
write.amelia(obj=a.out, file.stem = "a.out")


# Now, we just want to estimate a basic regression model with our imputed data
# But we have 5 datasets!
# Zelif will combine them for us:
z.out.imp <- zelig(trade ~ log(population) + log(gdp_pc) + infl 
                   + civlib, data = a.out$imputations, model = "ls") 
summary(z.out.imp)

### Diagnostics

# The missingness map shows us where our missingness is
missmap(a.out)

# We can plot our amelia output to compare empirical and imputed densities
plot(a.out)

# There are time-series, cross-sectional plots with imputed values
tscsPlot(a.out, var = "trade", cs = "Cameroon")

#Overimputation for a specific variable tests the imputation model by 
# imagining that each observation is missing and generating some imputations to check performance.   
overimpute(a.out, var = "trade")


#The disperse function starts the algorithm at some unlikely values to check 
# that amelia hasn't found a local rather a global maximum for the likelihood of the complete data.  
disperse(a.out, dims = 1, m = 5)





