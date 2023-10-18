
dt <- seq(0.05,0.15,0.02)
lambda = 1000 
deltaT = 0.05
deltaI = log(2)/1.39
beta<-5/(parmsTI["lambda"]/parmsTI["deltaT"]/parmsTI["deltaI"])

xstartTI.Manual <- c(Ta = lambda/deltaT,
                     I = 1)


for (t in seq(length(dt))){
  #the dynamics of the systme is store in the following matrix
  out.TI.manual <- rbind(c(t=0, Ta = xstartTI.Manual["Ta"], I = xstartTI.Manual["I"]))
  
  for(k in 1:1000){
    #(1) extract the current state of the system from the matrix
    t  <- out.TI.manual[nrow(out.TI.manual), "t"]
    Ta <- out.TI.manual[nrow(out.TI.manual), "Ta"]
    I  <- out.TI.manual[nrow(out.TI.manual), "I"]
    
    #Calculate the increments (i.e. the changes of the variables) that occur in the time-step dt
    # 
    dT <- (lambda *deltaT* Ta  - beta * Ta * I)* dt[t]
    dI <-  (beta * Ta * I - deltaI * I )* dt[t]
    
    #calculate the state of the system after this time step
    t_updated  = t + dt[t]
    Ta_updated = Ta + dT
    I_updated  = I + dI
    
    #add this to the matrix
    out.TI.manual <- rbind(out.TI.manual, c(t_updated, Ta_updated, I_updated))
  }
}


out.TI.manual<-as.data.frame(out.TI.manual)

plot(out.TI.manual)
