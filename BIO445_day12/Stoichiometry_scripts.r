############################################################
#############                                  #############
#############       PK/PD & Stoichiometry      #############
#############                                  #############
############################################################

############# Stoichiometry scripts


# Sample a trimer ---------------------------------------------------------

sample_trimer <- function(f_M) {

  # "M": mutant Env protein
  # "W": wild-type Env protein
  proteins <- c("M", "W") # possible proteins to sample from

  # randomly sample 3 proteins with replacement and prob=f_M
  trimer <- sample(proteins,
                   size = 3,
                   replace = TRUE,
                   prob = c(f_M, 1-f_M))

    return(trimer)

}
trimer_1 <- sample_trimer(0.5)
trimer_2 <- sample_trimer(0)

is_functional <- function(trimer){
  return (all(trimer %in% c("W")))
}

is_functional(trimer_1)
is_functional(trimer_2)

# Sample a virion ---------------------------------------------------------

sample_virion <- function(s, f_M) {

  virion_trimers <- replicate(n = s,
                              expr = sample_trimer(f_M),
                              simplify = FALSE)
  return(virion_trimers)

}

virion_1 <- sample_virion(4,0.1)


# number of functional trimer on virion
g_trimers <- function(virion){
  return (sum(sapply(virion, is_functional)))
}

g_trimers(virion_1)


# sample 10^4 virions with 13 trimers and f_M = 0.103
s <- 13
f_M <- 0.103
virions <- replicate(n = 10000,
                      expr = sample_virion(s,f_M),
                      simplify = FALSE)

# distribution of number of functional trimers
nr_func_trimers <- unlist(lapply(X = virions, FUN = g_trimers))

# histogram of functiounal primers
h <- hist(nr_func_trimers,
          col = "lightgray",
          border = "white",
          freq = FALSE,
          xaxt="n",
          main = "Distribution of the number of functional trimers",
          xlab = "Number of functional trimers",
          ylab = "Probability",
          breaks = seq(0,s+1)
          )
axis(side = 1, at=h$mids[1:s], labels = h$breaks[2:(s+1)])
# theoretical distribution
p <- (1-f_M)^3
pdf_binom <- dbinom(x = seq(1,s), size = s, prob = p)
lines(h$mids[1:s],pdf_binom, lty="dotted")
points(h$mids[1:s], pdf_binom, pch=19,cex=1.2)


# infectivity
infectivity <- function(virions, TT){
  nr_func_trimers <- unlist(lapply(X = virions, FUN = g_trimers))
  
  inf_virions <- which(nr_func_trimers > TT)
  
  prop_inf <- length(inf_virions)/length(nr_func_trimers)
  
  return (prop_inf)
}

# simulate different TT

prop_inf <- sapply(seq(1,14),function(x){
  return(infectivity(virions,x))
})

plot(prop_inf)


# for wildtype
virions_w <- replicate(n = 10000,
                     expr = sample_virion(s,0),
                     simplify = FALSE)
prop_inf_w <- sapply(seq(1,14),function(x){
  return(infectivity(virions_w,x))
})

plot(prop_inf_w)



# import empricial trimer numbers
library(readr)
trimer_numbers <- read_csv("Uni/Computational/HS22/BIO445/BIO445_day12/trimer_numbers.csv")

hist(trimer_numbers$mean_eta)

mean_trimer <- mean(trimer_numbers$mean_eta)

# experimental relative infectivity
RI_data <- read_csv("Uni/Computational/HS22/BIO445/BIO445_day12/RI_data.csv")
library(ggplot2)

ggplot(RI_data, aes(x = as.factor(f.M), y = RI, color = mutation))+
  geom_point() +
  facet_wrap(~virus) +
  geom_smooth()


# create a dataset for Cap88
library(dplyr)
cap88_data_wt <- RI_data %>%
  filter(virus == "Cap88" & mutation == "V513E" & env == "wt")

cap88_data_mut <- RI_data %>%
  filter(virus == "Cap88" & mutation == "V513E" & env == "V1V2")


ggplot(cap88_data_mut, aes(x = f.M, y = RI))+
  geom_point()

ggplot(cap88_data_wt, aes(x = f.M, y = RI))+
  geom_point()


# Stoichiometry estimator -------------------------------------------------

estimate_T <- function(data, trimer_number_sample) {

  ## Make sure that the data format is data.frame
  data <- as.data.frame(data)

  ## Fit trimer number distribution
  pdf_eta <- function(trimer_number_sample) {

    # trimer_number_sample = vector containing trimer numbers

    mu_s <- mean(trimer_number_sample)
    var_s <- 49/14*mu_s
    max_s <- 100

    #### definition of the trimer number distribution
    # correction for infinity
    eta_muv.b<-function(mu,v,s.max){
      # mu:		mean of trimer number distribution
      # v:		variance of trimer number distribution
      # s.max:	maximal number of trimers on virion surface (minimal number is 0)

      # calculation of parameters of B-distribtion

      mu<-mu/s.max
      v<-v/s.max^2

      p<-(mu^2-mu^3-mu*v)/v
      q<-(mu-2*mu^2+mu^3-v+v*mu)/v
      dos<-dbeta((0:s.max)/s.max,p,q)

      if(dos[1]==Inf){dos[1]<-ceiling(dos[2])
      }
      if(dos[s.max+1]==Inf){dos[s.max+1]<-ceiling(dos[s.max])
      }

      out<-dos/sum(dos)

    }

    distr <- eta_muv.b(mu_s, var_s, max_s)
    names(distr) <- seq(0, max_s, by=1)

    return(distr)

  }

  eta <- pdf_eta(trimer_number_sample)

  # RI function from the basic model
  entryRI.basic<-function(eta,f.M,TT){

    s.max<-length(eta)-1
    p3<-(1-f.M)^3

    if(f.M==0 | f.M==1){out<-1-f.M}

    else{
      out<-sum(eta[(TT:s.max)+1]*sapply(TT:s.max,function(s) sum(dbinom(TT:s,s,p3))))/sum(eta[(TT:s.max)+1])
    }

    out

  }

  # vectorization
  entryRI.basic.v<-function(eta,F.M,TT)sapply(F.M, entryRI.basic,eta=eta,TT=TT)

  # arcsin sqrt - transformation
  trafas<-function(x){
    if(!is.na(x)){
      if(x<0){x<-x}
      else{x<-asin(sqrt(x))}
    }
    x
  }
  arcsinsqrt<-function(x) sapply(x,trafas)

  # rss. matrix
  rss.matrix.entry.basic<-function(eta,data){
    T.max<-length(eta)-sum(eta==0)
    #	print(paste("T.max",T.max))
    out<-matrix(NA,nrow=T.max,ncol=2,dimnames=list(rep("",T.max),c("TT","rss")))
    for(TT in 1:T.max){

      #print(paste("T=",TT,sep=""))

      # no transformation
      # out[TT,]<-c(TT,sum(((data[,"RI"])-(entryRI.basic.v(eta,data[,"f.M"],TT)))^2))
      # arcsin sqrt- transformation
      out[TT,]<-c(TT,sum((arcsinsqrt(data[,"RI"])-arcsinsqrt(entryRI.basic.v(eta,data[,"f.M"],TT)))^2))
      # logit transformation
      # out[TT,]<-c(TT,sum((logit(data[,"RI"])-logit(entryRI.basic.v(eta,data[,"f.M"],TT)))^2))
    }
    out
  }

  # estimator
  T.estimator.entry.basic<-function(eta,data){
    W<-rss.matrix.entry.basic(eta,data)
    array(c(W[which(W[,"rss"]==min(W[,"rss"],na.rm=TRUE)),c("TT","rss")]),dim=c(1,2),dimnames=list("",c("TT","rss")))
  }

  T_est <- T.estimator.entry.basic(eta, data)
  return(T_est)
}

# prepare trimer numbers

cap88_trimer_mut <- trimer_numbers %>%
  filter(virus == "Cap88" & env == "V1V2")

cap88_trimer_wt <- trimer_numbers %>%
  filter(virus == "Cap88" & env == "wt")
  
estimate_T(cap88_data_mut,cap88_trimer_mut$mean_eta)
estimate_T(cap88_data_wt,cap88_trimer_wt$mean_eta)

# slide
virion_1 <- sample_virion(4,0.1)


g_trimers(virion_1)


# sample 10^4 virions with 13 trimers and f_M = 0.103

x = unique(RI_data$virus)[6]

create <- function(x){
  trimer_mut <- trimer_numbers %>%
    filter(virus == x & env == "V1V2")
  nr_trimer_mut <- mean(trimer_mut$mean_eta)
    
  trimer_wt <- trimer_numbers %>%
    filter(virus == x & env == "wt")
  nr_trimer_wt <- mean(trimer_wt$mean_eta)
  
  f_M <- 0
  
  # compute virions with corresponding number of trimers
  nr <- 1000
  virions_wt <- rpois(nr, nr_trimer_wt)
  
  virions_mut <- rpois(nr, nr_trimer_mut)
  
  
  # estimate TT
  RI_wt <- RI_data %>%
    filter(virus == x & mutation == "V513E" & env == "wt")
  
  RI_mut <- RI_data %>%
    filter(virus == x & mutation == "V513E" & env == "V1V2")
  
  estimated_TT_mut <- estimate_T(RI_mut,nr_trimer_mut)
  estimated_TT_wt <- estimate_T(RI_wt,nr_trimer_wt)
  
  
  mut_infectivity <- sum(virions_mut >= estimated_TT_mut[1]) / nr
  wt_infectivity <- sum(virions_wt >= estimated_TT_wt[1]) / nr

  df <- data.frame(nr_trimers = c(virions_mut,virions_wt), env = c(rep("mut",length(virions_mut)),rep("wt",length(virions_wt))), estimated_TT = c(rep(estimated_TT_mut[1]),rep(estimated_TT_wt[1])))
  
  ggplot(df%>%filter(env == "mut"),aes(x = nr_trimers))+
    geom_density() +
    ggtitle(paste0("strain: ",x,", estimated TT:",estimated_TT_mut, ", infectivity: ",mut_infectivity)) + 
    xlim(0,max(df$nr_trimers))+
    geom_vline(xintercept = estimated_TT_mut[1], color = "red")
  
  ggplot(df%>%filter(env == "wt"),aes(x = nr_trimers))+
    geom_density() +
    ggtitle(paste0("strain: ",x,", estimated TT:",estimated_TT_wt, ", infectivity: ",wt_infectivity)) + 
    xlim(0,max(df$nr_trimers))+
    geom_vline(xintercept = estimated_TT_wt[1], color = "red")
  
}


sapply(unique(RI_data$virus), create)



