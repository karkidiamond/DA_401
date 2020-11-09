library(tidyverse)
library(dplyr)
library(purrr)
library(stringr)
library(ggfortify)
library(caret)
library(plotly)
library(dplyr)
library(mice)
library(caret)
library(randomForest)
library(mlbench)
library(e1071)
library(ggplot2)



morgan<-read.csv('final_morgan.csv')
#removing Zero variance column
removeZeroVar2 <- function(df){
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}

morgan<-removeZeroVar2(morgan)#eliminated 53 zero variance columns

fit_data<-morgan[-c(1:39,42,44)]

#fit_data <- fit_data[!is.na(fit_data$ORAL_RAT_LD50_MOL.KG_TEST_PRED), ]
sub_fit_data<-fit_data[1:2,]


#cross-validation for splitting
#5-fold cross validation 
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 2,
  savePredictions = TRUE,
  classProbs = FALSE,
  verboseIter = TRUE,
  preProcOptions=list(thresh=0.95,na.remove=TRUE,verbose=TRUE))


trees<-c(1000,1500,2000)
rf_model_list<-list()#initate list of models
set.seed(12345)
for (n in trees){
  rf_model <- train(xlogp ~ .,
                    data=fit_data,
                    method="rf",
                    na.action = na.exclude,
                    preProcess=c("center","scale"),#scale the data
                    tuneLength=12,
                    metric="RMSE",
                    trControl=ctrl,
                    ntree = n)#different trees
  key <-toString(n)
  rf_model_list[[key]] <-rf_model
  
}

#run time 9 hrs 34 min 5s

rf_best<-as.data.frame(rf_model_list[[2]]$results)

#rf_model with 1000 tress produced lowest RMSE=0.9217127,Rsquared=0.9141306 with mtry=830
#rf_model with 1500 tress produced lowest RMSE=0.9055821,Rsquared=0.9123824 with mtry=1623
#rf_model with 2000 tress produced lowest RMSE=0.9071887,Rsquared=0.9123356 with mtry=1623




fit_data<-removeZeroVar2(fit_data)#eliminate zero variance columns
#install.packages("logisticPCA")
#library(logisticPCA)
#library(rARPACK)
#logpca_cv = cv.lpca(morgan[1:3,64:1034], ks = 1:30, ms = 1:10)
# k <- 5; m <- 10
# logpca_model = logisticPCA(binary_df, k = k, m = m)
# 
# logpca_features <- predict(logpca_model, binary_df, type = "PCs") 
# colnames(logpca_features) <- paste0("LPC", 1:k)

#Principle component to reduce dimensions/feature selection
pca_fit_data<-fit_data[,-c(12:15)]
fit.pca = prcomp(pca_fit_data[,-c(18)],scale. = TRUE) #perform pca on the data
pca_var = fit.pca$sdev^2
prop_varex <- pca_var/sum(pca_var) #Proportion of variance
which.max(cumsum(prop_varex)>.95)  #pc's needed for 95% of the variance

pc_data = data.frame(fit.pca$x[,1:303]) #keep only the 303 pc that explian 95% variance of the data
pc_data['xlogp']=fit_data$xlogp#adding xlogp to the pc data

#model with pc data of morgan fingerprints
pc_rf_model <- train(xlogp ~ .,
                     data=pc_data,
                     method="rf",
                     na.action = na.exclude,
                     tuneLength=12,
                     metric="RMSE",
                     trControl=ctrl,
                     ntree = 1500)

pc_rf_model

#visulaizing mtry vs RMSE
#Basic line plot with points
ggplot(data=pc_rf_model$results, aes(x=mtry, y=RMSE, group=1)) +
  geom_line()+
  geom_point()
#Change the line type

pc_rf_model$pred

#rm(basic_fit)
#Model with only basic molecular descriptors form pubchem, no morgan fingerprints
basic_fit<-fit_data[,-c(12:15,23:993)]
basic_rf_model <- train(xlogp ~ .,
                        data=basic_fit,
                        method="rf",
                        na.action = na.exclude,
                        mtry=1623,
                        metric="RMSE",
                        trControl=ctrl,
                        ntree = 1500)
basic_rf_model
