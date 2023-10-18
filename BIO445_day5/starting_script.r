##### STARTING SCRIPT - ECOLOGICAL NETWORKS

# SETUP #
rm(list=ls())

#libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rjson, data.table, bipartite, igraph, Matrix)





##### 1. Network Science -----------------------------------------------

# A) ----------
download_and_write_csv_single_network <- function(network_name){
  networkName <- network_name
  speciesName <- "yes"
  url <- paste("http://www.web-of-life.es/download/", networkName, "_", speciesName, ".csv", sep = "")
  data <- fread(url)
  write.csv(data, file = paste(trimws(networkName), ".csv", sep = ""))
  return(data)
}


network <- download_and_write_csv_single_network("M_PL_052")
names <- network$V1
network <- as.matrix(network[,-1])

rownames(network) <- names

# B) ----------

library(ComplexHeatmap)

Heatmap(network)

plotweb(network)



# C) ----------

data <- network
n_row <- nrow(network)
n_col <- ncol(network)
full_adjm <- function(data, n_row,n_col) {
  full_adjm <- matrix(0, n_row + n_col, n_row + n_col)
  full_adjm[(n_row + 1):(n_row + n_col), 1:n_row] = t(as.matrix(data > 0))
  full_adjm[1:n_row, (n_row + 1):(n_row + n_col)] = as.matrix(data > 0)
  return(full_adjm)
}

adjm <- full_adjm(data, n_row, n_col)
net <- graph_from_adjacency_matrix(adjm, mode = "undirected")
colrs <- c("tomato", "gold")
V(net)$color <- colrs[as.numeric(V(net) > n_row) + 1]
V(net)$label <- append(row.names(data), colnames(data))
plot(net)
plot(net, layout = layout_in_circle)



# D) ----------
library(ggplot2)
library(dplyr)
library(tidyverse)

rowdegree <- rowSums(network)
coldegree <- colSums(network)

df <- data_frame(species = c(names(rowdegree),names(coldegree)), degree = c(rowdegree, coldegree)) %>%
  arrange(degree)
  


ggplot(df, aes(x = as.factor(species), y = degree)) + 
  geom_point()



# E) ----------
# Hint: Density is the ratio of observed interactions to all possible interactions in a network



density <- function(matrix){
  # calculate connentedness here and store it as the variable result
  result <- sum(rowSums(matrix))/ (nrow(matrix)*ncol(matrix))
  return(result)
}

densi <- density(network)



# F) ----------
# Hint: write a function that returns the nestedness

nestedness_binmatnest2 <- function(matrix) {
  # The following function calculates the nestedness temperature, which ranges between 0 (fully nested) and 100 (not nested)
  value <- as.numeric(nested(matrix, "binmatnest"))
 
  # modify value here as indicated in the exercise sheet
  value <- (100-value)/100
  
  return(value)
}


nest <- nestedness_binmatnest2(network)




# G) ----------
# Hint: write a function that takes a matrix as input. 
#       First generate an empty matrix of the same size.
#       Next read in Prof. Bascompte´s paper, how the empty matrix will be filled in, depending on the input matrix

# Hint: random number generators work as follows:

nestedness_list <- c()
density_list <- c()
for (i in 1:1000){
  sample <- rbinom(n_col*n_row, 1, prob = densi)  # returns 5 values that are integers <= 1, with the probability per trial of 0.6
  mat <- matrix(sample,nrow = n_row,ncol = n_col)
  if (any(rowSums(mat) >= 1) & any(colSums(mat) >= 1)){
    nested <- nestedness_binmatnest2(mat)
    nestedness_list <- c(nestedness_list,nested)
    dens <- density(mat)
    density_list <- c(density_list,dens)
  }

}

hist(nestedness_list)
hist(density_list)




results <- sum /length(1:1000)



# H) ----------
# Hint: you have a function to calculate nestedness, and you have a function that generates you a random nullmodel. 



# I) ----------





##### 2. Scale-free networks -----------------------------------------------

# B) ---------- 
# create a network with preferential attachement
edgelist <- c(1, 2) # list of the edges
degreelist <- c(1, 1) # list of the degrees of each node

# Hint: use a loop. You already have two connected nodes, you will start with the third
#       sample.int(number_of_existing_nodes, number_of_edges, replace = F, prob = degreelist) should help you. start off with allowing only 1 edge per new node

A <- make_graph(edgelist, directed = FALSE)
V(A)$label <- ""
V(A)$size <- degreelist + 2
plot(A)
