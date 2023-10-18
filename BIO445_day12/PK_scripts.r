############################################################
#############                                  #############
#############       PK/PD & Stoichiometry      #############
#############                                  #############
############################################################

############# Pharmacokinetics scripts


# Lean body mass index ----------------------------------------------------

LBM <- function(gender, height, weight) {
  
  #### Arguments:
  # - gender: string indicating gender
  # - height: in centimeters
  # - weight: in kilograms
  
  if (gender == "M") {
    lbm <- 1.1*weight-128*(weight/height)^2
  }
  else if (gender == "F") {
    lbm <- 1.07*weight-148*(weight/height)^2
  } else {
    stop("Gender should be either 'M' for males or 'F' for females!")
  }
  return(lbm)
}


# Schnider model ----------------------------------------------------------
Schnider_model <- function(gender, height, weight, age) {
  
  # LBM
  lbm <- LBM(gender, height, weight)
  
  # volumes (L)
  V_1 <- 4.27
  V_2 <- 18.9-0.391*(age-52)
  V_3 <- 238
  
  # clearances (L/min)
  CL <- 1.89+0.0456*(weight-77)-0.0681*(lbm-59)+0.0264*(height-177)
  CL_2 <- 1.29-0.024*(age-52)
  CL_3 <- 0.836
  
  pars <- list("V_1"=V_1, "V_2"=V_2, "V_3"=V_3,
            "CL"=CL, "CL_2"=CL_2, "CL_3"=CL_3)
  return(pars)
}


# Single compartment - loading dose ---------------------------------------
c1_loading <- function(T_max, dt, V_1, CL, r, D_0, 
                       T_stop = T_max+dt) {
  
  #### T_max: maximum observed time point
  #### dt: time-stp for numerical solution of the ODE
  #### V_1: volume of distribution 1
  #### CL: clearance
  #### r: infusion rate
  #### D_0: loading dose
  #### T_stop: time of the infusion stop
  
  # model parameters: V_1, CL and r
  params <- c("V_1"=V_1,
              "CL"=CL,
              "r"=r)
  
  # initial concentrantions - 1 (+ 1 compartment)
  c_init <- c("c_1"=D_0/V_1, # loading dose concentration
              "c_1_memory"=D_0/V_1) # serves as a memory of c_1
  
  # c_1_memory is needed to memorize the concentration at the infusion stop
  # when the infusion is stopped the derivative of the c_1_memory will be changed to 0
  # so that the c_1_memory concentration stays the same as at the time point
  # of stopping the infusion
  # To find the context sensitive half-life we will aim to detect the event when
  # the c_1_memory concentration is twice a high as c_1
  
  # time points
  times <- seq(from=0, to=T_max, by=dt)
  
  # specify the derivative functions
  dc <- function(time, conc, params) {
    
    # what is the infusion rate at time:
    rate <- params["r"]
    if (time >= T_stop) { # after the infusion stop
      rate <- 0
    }
    
    # dc_1
    dc_1 <- 1/params["V_1"]*(rate-params["CL"]*conc["c_1"])
    
    # memory of c_1
    if (time < T_stop) { # before the infusion stop -> same as for c_1
      dc_1_memory <- dc_1
    }
    else { # after the infusion stop
      dc_1_memory <- 0 # stops changing so that it stays the the level of infusion stop
    }
    
    derivatives <- c("c_1"=dc_1,
                     "c_1_memory"=dc_1_memory)
    
    return(list(derivatives))
  }
  
  # concentration hit function: 
  # when the specified concentration (c_target) is reached
  # in our case c_target = 0.5*conc["c_1_memory"] since after the infusion stop
  # c_1_memory concentrations will remain unchanged
  hit_c <- function(time, conc, params) {
    dif <- conc["c_1"] - 0.5*conc["c_1_memory"]
    return(dif)
  }
  
  # are the concentrations modified upon event?
  action <- function(time, conc, params) {
    ### if for instance c_1 is increased by fixed value c_increase
    # conc["c_1"] <- conc["c_1"] + c_increase
    # return(conc)
    
    ### no modification (we just need the time of the event)
    return(conc)
  }
  
  # numerical solution
  c_out <- ode(y = c_init,
               times = times,
               func = dc,
               parms = params,
               events = list(func=action, # how are the concentrations affected upon the event?
                             root=TRUE), # the event is of a type: when a certain level is crossed
               rootfun = hit_c) # function measuring how far away from crossing the level are the current concentrations
  
  # time to event
  t_t_e <- attributes(c_out)$troot # how to extract time of the event
  attr(c_out, "time_to_event") <- t_t_e
  
  return(c_out)
  
}

