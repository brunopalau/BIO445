
# Problem 2 ----
## 2.1

strategy1 <- function(ownpastactions, partnerpastactions) { return(TRUE) } # always TRUE

strategy2 <- function(ownpastactions, partnerpastactions) { return(runif(1) > 0.75) } # 1/4 TRUE

strategy3 <- function(ownpastactions, partnerpastactions) { # 1/3 TRUE
    n <- length(ownpastactions)
    
    if (n %% 3 == 0) {
        return(TRUE) 
    } else { 
        return(FALSE) 
    } 
}

stategyBruno <- function(ownpastactions, partnerpastactions){
  short_term <- 100
  # compute long term probability of pascal to say TRUE
  long_term <- sum(partnerpastactions)/length(partnerpastactions)
  # compute short term probability of pascal to say TRUE
  if (length(partnerpastactions) > 5){
    short_term <- sum(partnerpastactions[1:5])/5
  }
  
  # otherwise choose average of both
  if (short_term != 100){
    # if pascal only says false in last steps punish him till he says TRUE
    if (short_term == 0){return (FALSE)}
    # if he always says TURE, very often say TRUE
    if (short_term == 1){return (rnorm(1,0) > -2)}
  return(runif(1) < ((short_term + long_term)/2))
  }
  return (runif(1) < long_term)
}

strategyPascal1 <- function(ownpastactions, partnerpastactions){ # tit for tat
  n <- length(ownpastactions)
  i <- c(1:(n - 1))
  
  if (n == 0){return (T)}
  
  else {return (partnerpastactions[-i])}
  
}

strategyModel <- function(ownpastactions,partnerpastactions){
  
  # create feature data frame
  df <- data.frame(correlation = c(), total_true = c(), after_false = c(), after_true = c(), response = c())
  
  if (length(ownpastactions) < 10){
    start <- 1
  }
  else{start <- length(partnerpastactions)-10}
  
  new_entry <- c(cor(ownpastactions,partnerpastactions), sum(partnerpastactions[start:length(partnerpastactions)]),
                 sum(partnerpastactions[which(ownpastactions[1:length(ownpastactions)-1] == F)+1]),
                 sum(partnerpastactions[which(ownpastactions[1:length(ownpastactions)-1] == T)+1]),
                 partnerpastactions[length(partnerpastactions)])
  
  # add new entry to dataframe
  df <- rbind(df,new_entry)
  colnames(df) <- c("correlation", "total_true", "after_false", "after_true", "response")
  
  model <-  glm(response ~.,family=binomial(link='logit'),data=df)
  
  predict(model, new_entry)
  
}



## 2.2

payout <- function(action_p1, action_p2, payoutvector = c(5, 3, 1, 0)) {
    # best (5) > good (3) > bad (1) > worst (0)
    if (action_p1 & action_p2) return(c(payoutvector[2], payoutvector[2])) 
    if (action_p1 & ! action_p2) return(c(payoutvector[4], payoutvector[1])) 
    if (! action_p1 & action_p2) return(c(payoutvector[1], payoutvector[4])) 
    if (! action_p1 & ! action_p2) return(c(payoutvector[3], payoutvector[3]))
}

simulate_2P_PD <- function(strategy_p1, strategy_p2, N, payoutvector = c(5, 3, 1, 0)) { 
    total_payout <- c(0, 0)
    past_actions_p1 <- vector()
    past_actions_p2 <- vector()
    for (i in 1:N) { 
      
      decision_p1 <- strategy_p1(past_actions_p1,past_actions_p2)
      decision_p2 <- strategy_p2(past_actions_p2,past_actions_p1)
      out <- payout(decision_p1,decision_p2)
      
      past_actions_p1 <- c(past_actions_p1, decision_p1)
      past_actions_p2 <- c(past_actions_p2, decision_p2)
      total_payout <- total_payout + out
    }
    return(total_payout) 
}


simulate_2P_PD(strategy1,strategy2,2)

simulate_2P_PD(stategyBruno,strategyPascal1,50)

## 2.3

tournament <- function(list_of_strategies = list(strategy1, strategy2, strategy3), 
                       N = 50, Nsim = 100, payoutvector = c(5, 3, 1, 0), 
                       PLAY_YOURSELF = FALSE) {
    
    n <- length(list_of_strategies) 
    total_payouts = rep(0, n) 
    av_payout_mat = matrix(0, n, n)
    
    for (i in 1:(n - 1 + PLAY_YOURSELF)) {
        for (j in (i + 1 - PLAY_YOURSELF):n) {
            average_payout <- ?? #Remember to average 100 repetitions 
                av_payout_mat[i,j] <- ??
                    av_payout_mat[j,i] <- ??
                        total_payouts[i] <- total_payouts[i] + average_payout[1] / 
                            (n - 1 + PLAY_YOURSELF) 
                        total_payouts[j] <- total_payouts[j] + average_payout[2] / 
                            (n - 1 + PLAY_YOURSELF)
        } 
    }
    return(av_payout_mat) 
}

heatmap(av_payout_mat, scale = 'none', revC = TRUE, Rowv = NA, Colv = NA, margins = c(10, 10))

