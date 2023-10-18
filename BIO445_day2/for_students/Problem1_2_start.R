library(deSolve)
library(ggplot2)

TI.model <- function(t, x, parms) {
  Ta <- x[1] # pop 1
  I  <- x[2] # pop 2
  with(as.list(parms),{
    dT <- lambda- deltaT * Ta - beta * Ta * I
    dI <- beta * Ta * I - deltaI * I 
    res <- c(dT, dI)
    list(res)
  })
}

### by varying beta
lambda <- seq(from=500,to=1500,by=50)
deltaT <- seq(from=0.01,0.1,0.01)
deltaI <- seq(from=log(2)/1.39 - 0.2,log(2)/1.39 + 0.2,0.02)
beta <- seq(4,6,0.2)


out_df <- data_frame()
for (i in c(1:10)){
  parmsTI  <- c(lambda = 1000, 
                deltaT = 0.05,
                deltaI = deltaI[i],
                beta = NA)
  
  parmsTI["beta"] <- 5/(parmsTI["lambda"]/parmsTI["deltaT"]/parmsTI["deltaI"])
  
  xstartTI <- c(Ta = as.numeric(parmsTI["lambda"]/parmsTI["deltaT"]), # why is Ta dependent of lambda and deltaT
                I = 1)
  times  <- seq(0, 100, length=500)
  
  #here rk4 is a function provided by R to solve differential equations which you do not have to understand in further detail
  out.TI <- as.data.frame(rk4(xstartTI, times, TI.model, parmsTI))
  out.TI$id <- rep(deltaI[i],nrow(out.TI))
  
  out_df <- rbind(out_df,out.TI)
}
  
ggplot(out_df,aes(x = time, y= Ta, col = as.factor(id)))+
  geom_point()+
  geom_point(aes(y = I, color = as.factor(id)))
  
