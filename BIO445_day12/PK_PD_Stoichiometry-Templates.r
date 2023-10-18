############################################################
#############                                  #############
#############       PK/PD & Stoichiometry      #############
#############                                  #############
############################################################

############# Templates

######################################################################
### Problem 1
######################################################################

### (3.b)
c3_infusion <- function(T_max, dt, r,
                        V_1, CL, V_2, CL_2, V_3, CL_3) {
  
  #### T_max: maximum observed time point
  #### dt: time-stp for numerical solution of the ODE
  #### V_1, V_2, V_3: volumes of distribution 1-3
  #### CL, CL_2, CL_3: clearances
  #### r: infusion rate
  
  # model parameters: 
  params <- c("V_1" = V_1,
              "CL" = ???,
              ???)
  
  # initial concentrantions (at time = 0)
  c_init <- c("c_1" = ???, # blood compartment
              "c_2" = ???, # brain compartment
              "c_3" = ???) # fat tissue compartment
  
  # time grid
  times <- seq(from=???, to=???, by=???)
  
  # instantaneous concentrations change
  # i.e. the derivatives describing the ODE system
  dc <- function(time, conc, params) {
    
    # infusion rate 
    rate <- params["r"]
    
    # dc_1/dt
    dc_1 <- 1/params["V_1"]*(rate+params["CL_2"]*conc["c_2"]+???)
    
    # dc_2/dt
    dc_2 <- ???
      
      # dc_3/dt
      ???
      
      derivatives <- c("c_1" = dc_1, #dc_1/dt
                       "c_2" = ???, #dc_2/dt
                       "c_3" = ???) #dc_3/dt
    
    return(list(derivatives))
  }
  
  # numerical solution
  c_out <- ode(y = ???, # initial values in each compartment
               times = ???, # time grid
               func = ???, #ODEs
               parms = ???) # model parameters
  
  return(c_out)
  
}

### (3.c)
c3_infusion_chl <- function(T_max, dt, r,
                        V_1, CL, V_2, CL_2, V_3, CL_3,
                        T_stop) {
  
  #### T_max: maximum observed time point
  #### dt: time-stp for numerical solution of the ODE
  #### V_1, V_2, V_3: volumes of distribution 1-3
  #### CL, CL_2, CL_3: clearances
  #### r: infusion rate
  
  # model parameters: 
  params <- ??? # check c3_infusion
  
  # initial concentrantions (at time = 0)
  c_init <- ??? # check c3_infusion
  ## add an artificial compartment that serves as a memory
  ## i.e., saves the concentration upon infusion cesation
  c_init["c_1_memory"] <- c_init["c_1"]
  
  # time grid
  times <- ??? # check c3_infusion
  
  # instantaneous concentrations change
  # i.e. the derivatives describing the ODE system
  dc <- function(time, conc, params) {
    
    # what is the infusion rate before the infusion stop?
    rate <- ???
    if (time ???) { # after the infusion stop
      rate <- ???
    }
    
    # dc_1/dt
    dc_1 <- ??? # check c3_loading
    
    # dc_2/dt
    dc_2 <- ??? # check c3_loading
      
    # dc_3/dt
    dc_3 <- ??? # check c3_loading
      
    # dc_1_memory/dt
    dc_1_memory <- ???
    if (time ???) { # after the infusion cesation
      dc_1_memory <- ???
    }
      
    derivatives <- ??? # checkk c3_loading
    # add the memory compartment
    derivatives["c_1_memory"] <- dc_1_memory
    
    return(list(derivatives))
  }
  
  # concentration hit function: 
  # when is the specified concentration (c_target) reached?
  hit_c <- function(time, conc, params) {
    ??? # check c1_loading function
  }
  
  # are the concentrations modified upon event?
  action <- function(time, conc, params) {
    ??? # check c1_loading function
  }
  
  # numerical solution
  c_out <- ode(y = ???, # initial values in each compartment
               times = ???, # time grid
               func = ???, #ODEs
               parms = ???, # model parameters
               events = list(func = ???, # how are the concentrations affected upon the event?
                             root = TRUE), # the event is of a type: when a certain level is crossed
               rootfun = ???) # definition of the event
  
  # time to event
  t_t_e <- attributes(c_out)$troot # how to extract time of the event
  attr(c_out, "time_to_event") <- t_t_e
  
  return(c_out)
  
}

######################################################################
### Problem 2
######################################################################

### (1.c)
E <- function(conc, E0, Emax, EC50, gamma) {
  tot_eff <- ???
  return(tot_eff)
}

######################################################################
### Problem 3
######################################################################

### (1.a.ii)
trimer1 <- ???
  trimer2 <- ???
  
  is_functional <- function(trimer) {
    
    ???
      
  }

is_functional(trimer1)
???
  
### (1.b.ii)
g_trimers <- function(virion) {
    
  # Which trimers are functional?
    
  # How many functional trimers are there?
}

### (2.a)
s <- ???
  f_M <- ???
  virions <- replicate(n = ???, # nr. of repetitions
                       expr = ???, # function call to be repeated n-times 
                       simplify = FALSE)
# Distribution of number of functional trimers
nr_func_trimers <- unlist(lapply(X = ???, FUN = ???))
## histogram
h <- hist(???, 
          col="lightgray", 
          border = "white",
          freq=FALSE,
          main="Distribution of the number of functional trimers",
          xaxt="n",
          xlab="Number of functional trimers",
          ylab="Probability",
          breaks=seq(0, s+1))
axis(side=1, at=h$mids[1:s], labels=h$breaks[2:(s+1)])
## theoretical distribution
p <- (1-f_M)^3
pdf_binom <- dbinom(x=seq(1, s), size=s, prob=p)
lines(h$mids[1:s], pdf_binom, lty="dotted")
points(h$mids[1:s], pdf_binom, pch=19, cex = 1.2)

### (2.b)
infectivity <- function(virions, TT) {
  
  # Number of functional trimers for each virion (see (a))
  nr_func_trimers <- unlist(lapply(X = ???, FUN = ???))
  
  # Which virions have enough functional trimers?
  inf_virions <- which(nr_func_trimers ???)
  
  # Proportion of "infectious" virions
  prop_inf <- ???
    
    return(prop_inf)
}
