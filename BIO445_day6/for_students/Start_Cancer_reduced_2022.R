# Modeling Cancer


library(deSolve)

# define the ODE 
Gompertz <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dN  <- r*N*log(K/N)
    
    return(list(c(dN)))
  })
}

pars  <- c(r=0.015, K=40 )   # parameters   
yini  <- c(N=0.001) # initial values
times <- seq(0, 2000, by = 1) # time in days

# call the ODE solver
out   <- as.data.frame(ode(yini, times, Gompertz, pars))

plot(out, type="l", ylab="Tumor volume")