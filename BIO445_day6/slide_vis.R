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
plot(out)


## c) Effect of therapy
# Now let's examine the effect of therapy (chemotherapy or immunotherapy) by adding a right hand side negative term



### c.1) with treatment
Gompertz_Treat1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity), then introduce treatment at day 1095 (3 years)
    # HINT: overwrite alpha before time is 1095
    if(Time < 1095){dN  <- r * N * log(K/N)}else{dN  <- r * N * log(K/N) - alpha*c*N}
    return(list(c(dN)))
  })
}
pars <- c(r = 0.015, K = 40, alpha = 0.7, c = 0.3)
yini <- c(N = 0.001)
times <- seq(0, 2000, by = 1)

# plot the tumor volume considering treatment:
# HINT:
# if you're using ggplot, have a look at facet_zoom from the ggforce package

out   <- as.data.frame(ode(yini, times, Gompertz_Treat1, pars))
plot(out)

### c.2)
# Elimination: how many days does it take for tumor elimination?
# calculate and visualize

for(i in 1000:1500){
  if((out$N[i]< 0.001)){
    print(out[i,])
    break
  }
}


### c.3) 
# drug concentration decreases with time
Gompertz_Treat2 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity),
    # then introduce treatment at day 1095 (3 years)
    # The drug concentration also has to descrease only
    # since 1095 (before it was not in the body..)
    if(Time < 1095){
      dN  <- r * N * log(K/N)
      dC <- 0
    }else{
      dN  <- r * N * log(K/N) - alpha*C*N
      dC <- -0.25*C}
    return(list(c(dN, dC)))
  })
}

pars <- c(r = 0.015, K = 40, alpha = 0.7)
yini <- c(N = 0.001, C = 0.3)
times <- seq(0, 2000, by = 1)

# Visualize Tumor dynamics considering the decreasing drug concentration
out   <- as.data.frame(ode(yini, times, Gompertz_Treat2, pars))
plot(out$N)


# 50% of the tumor volume after treatment initiation
for(i in 1000:1500){
  if((out$N[i]< (max(out$N)/2)+1) & (out$N[i]> (max(out$N)/2-1))){
    print(out[i,])
  }
}

# Elimination:
for(i in 1000:1500){
  if((out$N[i]< 0.001)){
    print(out[i,])
    break
  }
}












################################################################################

Gompertz_Treat3 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    # let the tumor develop (reach carrying capacity),
    # then introduce treatment at day 1095 (3 years)
    # The drug concentration also has to descrease only
    # since 1095 (before it was not in the body..)
    if(Time < 500){
      dN  <- r * N * log(K/N)
      dC <- 0
    }
    else{
      
      dN  <- r * N * log(K/N) - alpha*C*N
      dC <- -c*C
      if(round(Time) %in% c(600,700,800,900,1000,1100,1200)){
        dC <-0.3 - c*C
      }
      if(N < 1){dN = 0}
    }
    return(list(c(dN, dC)))
  })
}

params_list <- list()
for (a in round(seq(0,1,length.out=10))){
  for (c in round(seq(0.01,0.1,length.out=5))){
    params_list <- append(params_list,list(c(a,c)))
  }
}


out_all <- lapply(params_list,function(x){
  alpha = x[[1]]
  c = x[[2]]
  pars <- c(r = 0.015, K = 40, alpha = alpha, c=c)
  yini <- c(N = 0.001, C = 0.3)
  times <- seq(0, 1300, by = 1)
  
  out_first   <- as.data.frame(ode(yini, times, Gompertz_Treat3, pars))
  out <- out_first
  out$alpha_start <- rep(alpha,nrow(out_first))
  out$c_start <- rep(c, nrow(out_first))
  print(alpha)
  return (out)
})

df_all <- data.frame(time = c(),N = c(), C = c(), alpha_start = c(), c_start = c())
for (i in 1:length(out_all)){
  df <- out_all[[i]]
  df_all <- rbind(df_all, df)
}

df_all$c_start <- as.factor(df_all$c_start)

ggplot(df_all,
       aes(x = time,y = N, color = c_start))+
  geom_point()+
  facet_wrap(~alpha_start)




days_to_healthy <- lapply(out_all,function(x){
  times <- x %>%
    filter(N < 1 & time > 100)%>%
    select(time,c_start,alpha_start)
  
  if (nrow(times) == 0){
    return(data.frame(time = NA, c_start = x$c_start[1], alpha_start = x$alpha_start[1]))
  }
  return (times[1,])
})

days_df <- list()
for (i in 1:length(days_to_healthy)){
  df <- days_to_healthy[[i]]
  
  days_df <- rbind(days_df, df)
}




days_df <- days_df%>%
  mutate(treatments_till_healthy =as.factor(round((time-500)/100)+1))

ggplot(days_df, aes(x = alpha_start, y = c_start, color = treatments_till_healthy, size = 10))+
  geom_point()+
  scale_color_manual(breaks = c("1", "2", "3", "NA"),
                                 values=c("red", "orange", "yellow", "grey"))


















### c.4)
# increasing the drug effectivness
# increase the effectivness (alpha)
for(i in 1:10){
  pars <- c(r = 0.015, K = 40, alpha = 1/10*i, t1 = 1095, t2 =, t3 =) ## YOUR CODE HERE ##
  yini <- c(N = 0.001, C = 0.3)
  times <- seq(0, 2000, by = 1)
  
  out   <- as.data.frame(ode(yini, times, Gompertz_Treat2, pars))
  plot(out$N)
}


# Plot the tumor rebound and calculate time to 50% of the tumor volume after treatment initiation

## YOUR CODE HERE ##


# Elimination:
## YOUR CODE HERE ##
