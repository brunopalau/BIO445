# Modeling Cancer
# 15.12.21

# Remove current working environment
rm(list = ls())

if(!requireNamespace("pacman", quietly = T))   install.packages("pacman")
pacman::p_load("deSolve", "ggforce", "tidyverse")






# Exercise 1 - Estimate tumor growth rate ------------------------------------------------------



## a) Load the database into R and explore it
setwd("Uni/Computational/HS22/BIO445/BIO445_day6/for_students/")
db <- read.csv2("Nakamura_db_BT_start.csv")

# TASK: Explore and get familiarized with the dataset


## YOUR CODE HERE ##

library(dplyr)



## b) Calculations

db <- db %>%
  mutate(absolute_growth = Latest_tumor_vol - Initial_tumor_vol) %>%
  mutate(absolute_growth_rate_per_year = absolute_growth / (Follow.up.time/12))%>%
  mutate(relative.growth.rate = ((Latest_tumor_vol / Initial_tumor_vol)**(1/(Follow.up.time/12)) - 1)  * 100)


# relative growth year percent
db$relative.growth.rate =  ## YOUR CODE HERE ##
quantile(db$relative.growth.rate, probs = c(0.2, 0.8))
# FOR YOU TO COMPARE: if your relative growth-rate is correctly calculated, 
quantile(db$relative.growth.rate, probs = c(0.2, 0.8)) # should return 3.07455 25.21822 



## c) Visualization
# Make a plot of the age of the patient versus absolute growth rate.
# Is there an apparent correlation?

library(ggplot2)

ggplot(db, aes(x = Age, y= absolute_growth_rate_per_year))+
  geom_point() +
  geom_smooth(method = "lm")

ggplot(db, aes(x = Age, y= relative.growth.rate))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_smooth(method = "lm", formula = y ~ poly(x, 3))



ggplot(db, aes(x = Age, y = Initial_tumor_vol))+
  geom_point()+
  geom_smooth()


## YOUR CODE HERE ##



## d) Linear regression
# Construct a linear regression using lm(), with age and absolute growth rate and plot the regression fit.

## YOUR CODE HERE ##



## e) Calculate the tumor doubling
db$doubling.time =  (db$Follow.up.time/12) * (log(2)/log(db$Latest_tumor_vol/db$Initial_tumor_vol))
quantile(db$doubling.time, probs = c(0.2, 0.8))
# FOR YOU TO COMPARE: if your relative growth-rate is correctly calculated, 
# quantile(db$relative.growth.rate, probs = c(0.2, 0.8)) should return 3.082191 22.889495 



### e.1) 
# Do tumors that were larger at the first screen also tend to have higher doubling time?
# How can you test this hypothesis?

## YOUR CODE HERE ##

ggplot(db, aes(x = Initial_tumor_vol, y = doubling.time, ))+
  geom_point()+
  geom_smooth(method = "lm")

library(ggpubr)
ggplot(db, aes(x= Age, y = doubling.time, color = Calcification))+
  geom_point()+
  stat_compare_means(method = "t.test")


### e.2)
# Do calcified tumors have a higher doubling time?
# How can you test this hypothesis?
# Think of a way to visualize your conclusion.

## YOUR CODE HERE ##






# Exercise 3 - Modeling with ODEs ----------------------------------------



## a)
# define the ODE 
Gompertz <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN  <- r * N * log(K/N)
    return(list(c(dN)))
  })
}

pars  <- c(r = 0.015, K = 40)  # parameters   
yini  <- c(N = 0.001)          # initial values
times <- seq(0, 2000, by = 1)  # time in days



## b)
# calling the ODE solver
out   <- as.data.frame(ode(yini, times, Gompertz, pars))

# plot the tumor volume over time

## YOUR CODE HERE ##



## c) Effect of therapy
# Now let's examine the effect of therapy (chemotherapy or immunotherapy) by adding a right hand side negative term



### c.1) with treatment
Gompertz_Treat1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity), then introduce treatment at day 1095 (3 years)
    # HINT: overwrite alpha before time is 1095
    
    ## YOUR CODE HERE ##
    
    return(list(c(dN)))
  })
}
pars <- c(r = 0.015, K = 40, alpha = 0.7, c = 0.3)
yini <- c(N = 0.001)
times <- seq(0, 2000, by = 1)

# plot the tumor volume considering treatment:
# HINT:
# if you're using ggplot, have a look at facet_zoom from the ggforce package

## YOUR CODE HERE ##


# 50% of the tumor volume after treatment initiation

## YOUR CODE HERE ##



### c.2)
# Elimination: how many days does it take for tumor elimination?
# calculate and visualize

## YOUR CODE HERE ##



### c.3) 
# drug concentration decreases with time
Gompertz_Treat2 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity),
    # then introduce treatment at day 1095 (3 years)
    # The drug concentration also has to descrease only
    # since 1095 (before it was not in the body..)
    
    ## YOUR CODE HERE ##
   
    return(list(c(dN, dC)))
  })
}

pars <- c(r = 0.015, K = 40, alpha = 0.7)
yini <- c(N = 0.001, C = 0.3)
times <- seq(0, 2000, by = 1)

# Visualize Tumor dynamics considering the decreasing drug concentration
## YOUR CODE HERE ##


# 50% of the tumor volume after treatment initiation
## YOUR CODE HERE ##


# Elimination:
## YOUR CODE HERE ##



### c.4)
# increasing the drug effectivness
# increase the effectivness (alpha)
pars <- c(r = 0.015, K = 40, alpha = ## YOUR CODE HERE ## ) 
yini <- c(N = 0.001, C = 0.3)
times <- seq(0, 2000, by = 1)

# Plot the tumor rebound and calculate time to 50% of the tumor volume after treatment initiation

## YOUR CODE HERE ##


# Elimination:
## YOUR CODE HERE ##






# Exercise 4 - Resistance -------------------------------------------------------------

# Now let’s assume that 1% of the treated cells develop resistance to treatment.
# Modify the system accord- ingly (assume that c(t) is constant again).
Gompertz_Treat_Resistance1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity),
    # then introduce treatment at day 1095 (3 years)
    # in this the drug concentration also has to descrease only
    # since 1095 (before it was not in the body..)
    # let cells leaving the N compartment become resistant with rate m (CAVE: gompertz, not exponential)
    
    ## YOUR CODE HERE ##
    
    return(list(c(dN, dC, dR)))
  })
}



## a.) 
# One month after the first treatment, what would be the volume of the resistant cells?

## YOUR CODE HERE ##



## b.)
# Plot the dynamics of both cells types. When will the resistant cells predominate?
# Plot both strains

## YOUR CODE HERE ##
  


## c.*)
# Now assume that resistant cells can revert back and become susceptible again,
# with a certain rate ”q”, modify the system accordingly and explore.
Gompertz_Treat_Resistance2 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity),
    # then introduce treatment at day 1095 (3 years)
    # in this the drug concentration also has to descrease only
    # since 1095 (before it was not in the body..)
    # let cells leaving the N compartment become resistant with rate m (CAVE: gompertz, not exponential)
    # let resistant cells revert with rate q 
    
    ## YOUR CODE HERE ##
    
    return(list(c(dN, dC, dR)))
  })
}

# Plot both strains

## YOUR CODE HERE ##
