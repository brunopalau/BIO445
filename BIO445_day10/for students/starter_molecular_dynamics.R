#############################################
###### Molecular Dynamics
# prepared by Jessy Duran Ramirez
# 20.12.22
#############################################


# R preparations ---------------------------------------------------------------------------------
# Remove current working environment, load packages and set your working directory 
rm(list = ls())

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChemmineR")

library(ChemmineR)
library(caret)

setwd("Uni/Computational/HS22/BIO445/BIO445_day10/for students/")

# set seed for reproducibility
set.seed(123)

# Problem 1  Molecules in R -------------------------------------------------------------------------

## a) Load the information of tannin into R

tannin <- read.SDFset ("Tannin.sdf")
plot(tannin[[1]])



## b) Your favorite molecule!
## YOUR CODE HERE ##

tannin <- read.SDFset ("Structure2D_CID_117941915.sdf")
plot(tannin[[1]])




# Problem 2: Exploring potential candidates for the next drug  ---------------------------------------------

## a) Read in and explore the data set:
library(readr)
desc <-  read_csv("desc.csv")
library(ggplot2)

ggplot(desc, aes(x = Ki_nM,y = apol)) + geom_point()
  
## b) desc_unique: randomly delete duplicated IDs

shuffled_data= desc[sample(1:nrow(desc)), ]
##random select 1 experiment for each Monomer
descUnique1 <- shuffled_data[!(duplicated(desc$BindingDB.MonomerID)),]
  
## average multiple measurements
library(dplyr)
means <- aggregate(desc$Ki_nM, by=list(desc$BindingDB.MonomerID), mean)
means$BindingDB.MonomerID <- means$Group.1

descUnique2 <- merge(desc[!(duplicated(desc$BindingDB.MonomerID)),], means, by="BindingDB.MonomerID", all = TRUE) %>%
  select(-c(Group.1, Ki_nM))%>%
  rename(Ki_nM = x)


# continue with your subset:
## YOUR CODE HERE ##
  
  
## c) Look at correlation between the features
## YOUR CODE HERE ##


## d) linear regression

model <- lm(desc$Ki_nM ~ ., desc)

anova(model)
summary(model)


## e) training and test set
boolean <- rep(FALSE, nrow(descUnique2))
boolean[sample(1:nrow(descUnique2),round(0.7*nrow(descUnique2)))] <- TRUE

training <- descUnique2[boolean,]
test <- descUnique2[boolean,]

## f) fit models
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
fit1 <- train(Ki_nM ~ ., data = training, method = "rpart", trControl = fitControl)
fit2 <- train(Ki_nM ~ ., data = training, method = "lm", trControl = fitControl)
fit3 <- train(Ki_nM ~ ., data = training, method = "nnet", trControl = fitControl)

fit1$results

varImp(fit1)

summary(fit1)

## i) test set performance
pred1 <- predict(fit1,test)
pred2 <- predict(fit2,test)
pred3 <- predict(fit3,test)

postResample(pred1,test$Ki_nM)
postResample(pred2,test$Ki_nM)
postResample(pred3,test$Ki_nM)








