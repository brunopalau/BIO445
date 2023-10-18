############################################################
#############                                  #############
#############       PK/PD & Stoichiometry      #############
#############                                  #############
############################################################

############# Pharmacodynamics scripts


# PK blood-brain concentrations model -------------------------------------

conc_Blood_Brain <- function(T_max = 120, 
                             dt = 5/60, 
                             V_Bl = 5, V_Br = 1.2,
                             CL_Bl = 0.1, CL_Br = 0.3,
                             D0 = 10, r = 0.15,
                             T_stop = 75) {
  # T_max: the last observed time
  # dt: time-step
  # V_Bl: volume of distribution of the blood compartment
  # V_Br: volume of distribution of the brain
  # CL_Bl: clearance of the blood compartment (elimination)
  # CL_Br: clearance of the brain compartment 
  # D_0: initial bolus dose
  # r: infusion rate
  # T_stop: infusion duration
  
  # ODE
  dc <- function(t, conc, pars) {
    with(as.list(c(conc, pars)), {
      
      # dc1
      dc_bl <- 1/V_Bl*(r*(t<=T_stop) + CL_Br*c_Br - (CL_Bl+CL_Br)*c_Bl)
      
      # dc2
      dc_br <- 1/V_Br*(CL_Br*c_Bl - CL_Br*c_Br)
      
      list(c("c_Bl"=dc_bl, 
             "c_Br"=dc_br))
    })
  }
  
  # parameters
  params <- c("V_Bl"=V_Bl, "V_Br"=V_Br,
              "CL_Bl"=CL_Bl, "CL_Br"=CL_Br,
              "r" = r)
  
  # initial values
  init_c <- c("c_Bl"=D0/V_Bl,
              "c_Br"=0)
  
  # times
  times <- seq(0, T_max, by=dt)
  
  # solution
  c_out <- deSolve::ode(y = init_c, 
                        times = times, 
                        parms = params, 
                        func = dc)
  dimnames(c_out)[[2]] <- c("Time", "Blood", "Brain")
  return(c_out) 
}
