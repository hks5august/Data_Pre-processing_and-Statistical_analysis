#!/usr/bin/Rscript
#Usage:  Rcript all_regression.R input_data.csv 
args <- commandArgs(trailingOnly = TRUE)
library("rpart")
library("ipred")
library("survival")
library("survminer")
library("DStree")
library("party")
library("plyr")
library("RANN")
library("dplyr")
library("RSNNS")
library("kernlab")
# library("neuralnet")
# library("xgboost")
# library("brnn")
library("randomForestSRC")
library("randomForest")
library("mlbench") ### for MAE calculation
library("keras")
library("glmnet")
library("caret")
library(magrittr) # needs to be run every time you start R and want to use %>%
library(tidyverse)
library(Hmisc)
library(e1071)
library(hydroGOF)
library(randomForest)
#library(MASS)
library(caTools)

set.seed(7)
###### load file ######
##### first column contains sample_id, 2nd labels and 3rd onward features #######
mat1 = read.csv(args[1], header = T, sep =',')
#mat2 = mat1[2:ncol(mat1)]
######## Train/test split in 70:30 ratio ####
dt = sort(sample(nrow(mat1), nrow(mat1)*.8))
train_1<- mat2[dt,]
test_1<- mat2[-dt,]

train_id <- train_1[1]
valid_data1_id <- test_1[1]

########### create regression function for training data ##########################
reg_function <- as.formula(paste(names(train_1)[2], "~", paste(names(train_1)[3:ncol(train_1)], collapse=" + ")))

#f <- as.formula(paste(names(train_1)[1], "~", paste(names(train_1)[2:51], collapse=" + ")))
##################### 10-fold Cross Validation on Training data ############################## 
train.crossVal = trainControl(method = "cv", number = 10)

#################### Model Developement using svmRadial method at best parameters ################################
model_svr_rbf<- train(reg_function, data = train_1 , method = "svmRadial", trControl = train.crossVal)
##################### extract best parameters of the final model #######################
parameters_svr<-model_svr_rbf$bestTune
write.csv(parameters_svr,file = "Best_SVR_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
svr_rmse<-model_svr_rbf$results$RMSE
svr_r2<-model_svr_rbf$results$Rsquared
svr_mae<-model_svr_rbf$results$MAE
perf_svr<-cbind(svr_rmse,svr_r2,svr_mae)
write.csv(perf_svr,file = "perf_train_SVR",row.names=F, quote = F)
##################### save model #################################
save(model_svr_rbf, file = "SVR_model.rda")
##################### Prediction on validation dataset1 ########################
test1 = test_1[3:ncol(test_1)]
test_pred<-predict(model_svr_rbf, newdata = test1)

Prediction<-round(test_pred,2)
Prediction_result_valid2= cbind(valid_data2_id,OS_time_valid2)
write.csv(Prediction_result_valid2,file = "OS_time_Prediction_VH_vimp_valid2_svr",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_svr_rbf_valid2 = rmse(test_1[2],Prediction)
mae_svr_rbf_valid2 = MAE(test_1[2],Prediction)
r2_svr_rbf_valid2 = R2(test_1[2],Prediction)
cor_svr_rbf_valid2 = cor(test_1[2],Prediction)
test_results_valid2_svr<-cbind(cor_svr_rbf_valid2,rmse_svr_rbf_valid2,r2_svr_rbf_valid2,mae_svr_rbf_valid2)
write.csv(test_results_valid2_svr,file = "test_results_valid2_svr",row.names=F, quote = F)



############################## Random Forest Model ###################################
#################### Model Developement using svmRadial method at best parameters ################################
model_rf<- train(reg_function, data = train_1 , method = "rf", trControl = train.crossVal)

##################### save model #################################
save(model_rf, file = "vimp_rf_model.rda")

##################### extract best parameters of the final model #######################
parameters_rf<-model_rf$bestTune
write.csv(parameters_rf,file = "Best_RF_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
rf_rmse<-model_rf$results$RMSE
rf_r2<-model_rf$results$Rsquared
rf_mae<-model_rf$results$MAE
perf_rf<-cbind(rf_rmse,rf_r2,rf_mae)
write.csv(perf_rf,file = "perf_train_RF",row.names=F, quote = F)



##################### Prediction on test dataset ########################
test1_pred_rf<-predict(model_rf, newdata = test1)
####################### round upto 2 digit ##############################
OS_time_valid1_rf<-round(test1_pred_rf,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_rf= cbind(valid_data1_id,OS_time_valid1_rf)
write.csv(Prediction_result_valid1_rf,file = "OS_time_Prediction_VH_vimp_valid1_rf",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_rf_valid1 = rmse(test_1[2],OS_time_valid1_rf)
mae_rf_valid1 = MAE(test_1[2],OS_time_valid1_rf)
r2_rf_valid1 = R2(test_1[2],OS_time_valid1_rf)
cor_rf_valid1 = cor(test_1[2],OS_time_valid1_rf)

test_results_valid1_rf<-cbind(cor_rf_valid1,rmse_rf_valid1,r2_rf_valid1,mae_rf_valid1)
write.csv(test_results_valid1_rf,file = "test_results_valid1_rf",row.names=F, quote = F)
#######################################################################################################


############################## Linear Regression Model ###################################
#################### Model Developement using svmRadial method at best parameters ################################
model_LR<- train(reg_function, data = train1 , method = "lm", trControl = train.crossVal)

##################### save model #################################
save(model_LR, file = "vimp_LR_model.rda")
##################### extract best parameters of the final model #######################
parameters_LR<-model_LR$bestTune
write.csv(parameters_LR,file = "Best_LR_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
LR_rmse<-model_LR$results$RMSE
LR_r2<-model_LR$results$Rsquared
LR_mae<-model_LR$results$MAE
perf_LR<-cbind(LR_rmse,LR_r2,LR_mae)
write.csv(perf_LR,file = "perf_train_LR",row.names=F, quote = F)



##################### Prediction on test dataset ########################
test1_pred_LR<-predict(model_LR, newdata = valid_data1)
####################### round upto 2 digit ##############################
OS_time_valid1_LR<-round(test1_pred_LR,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_LR= cbind(valid_data1_id,OS_time_valid1_LR)
write.csv(Prediction_result_valid1_LR,file = "OS_time_Prediction_VH_vimp_valid1_LR",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_LR_valid1 = rmse(test_1[2],OS_time_valid1_LR)
mae_LR_valid1 = MAE(test_1[2],OS_time_valid1_LR)
r2_LR_valid1 = R2(test_1[2],OS_time_valid1_LR)
cor_LR_valid1 = cor(test_1[2],OS_time_valid1_LR)

test_results_valid1_LR<-cbind(cor_LR_valid1,rmse_LR_valid1,r2_LR_valid1,mae_LR_valid1)
write.csv(test_results_valid1_LR,file = "test_results_valid1_LR",row.names=F, quote = F)
#######################################################################################################



############################## BRNN Model ###################################
#################### Model Developement using svmRadial method at best parameters ################################
model_brnn<- train(reg_function, data = train1 , method = "brnn", trControl = train.crossVal)

##################### save model #################################
save(model_brnn, file = "vimp_brnn_model.rda")
##################### extract best parameters of the final model #######################
parameters_brnn<-model_brnn$bestTune
write.csv(parameters_brnn,file = "Best_brnn_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
brnn_rmse<-model_brnn$results$RMSE
brnn_r2<-model_brnn$results$Rsquared
brnn_mae<-model_brnn$results$MAE
perf_brnn<-cbind(brnn_rmse,brnn_r2,brnn_mae)
write.csv(perf_brnn,file = "perf_train_brnn",row.names=F, quote = F)


##################### Prediction on Validation dataset1 ########################
test1_pred_brnn<-predict(model_brnn, newdata = valid_data1)
####################### round upto 2 digit ##############################
OS_time_valid1_brnn<-round(test1_pred_brnn,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_brnn= cbind(valid_data1_id,OS_time_valid1_brnn)
write.csv(Prediction_result_valid1_brnn,file = "OS_time_Prediction_VH_vimp_valid1_brnn",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_brnn_valid1 = rmse(test_1[2],OS_time_valid1_brnn)
mae_brnn_valid1 = MAE(test_1[2],OS_time_valid1_brnn)
r2_brnn_valid1 = R2(test_1[2],OS_time_valid1_brnn)
cor_brnn_valid1 = cor(test_1[2],OS_time_valid1_brnn)

test_results_valid1_brnn<-cbind(cor_brnn_valid1,rmse_brnn_valid1,r2_brnn_valid1,mae_brnn_valid1)
write.csv(test_results_valid1_brnn,file = "test_results_valid1_brnn",row.names=F, quote = F)
#######################################################################################################


############################## KNN Model ###################################
#################### Model Developement using svmRadial method at best parameters ################################
model_knn<- train(reg_function, data = train1 , method = "knn", trControl = train.crossVal)

##################### save model #################################
save(model_knn, file = "vimp_knn_model.rda")
##################### extract best parameters of the final model #######################
parameters_knn<-model_knn$bestTune
write.csv(parameters_knn,file = "Best_knn_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
knn_rmse<-model_knn$results$RMSE
knn_r2<-model_knn$results$Rsquared
knn_mae<-model_knn$results$MAE
knn_LR<-cbind(knn_rmse,knn_r2,knn_mae)
write.csv(perf_knn,file = "perf_train_knn",row.names=F, quote = F)



##################### Prediction on Validation dataset1 ########################
test1_pred_knn<-predict(model_knn, newdata = valid_data1)
####################### round upto 2 digit ##############################
OS_time_valid1_knn<-round(test1_pred_knn,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_knn= cbind(valid_data1_id,OS_time_valid1_knn)
write.csv(Prediction_result_valid1_knn,file = "OS_time_Prediction_VH_vimp_valid1_knn",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_knn_valid1 = rmse(test_1[2],OS_time_valid1_knn)
mae_knn_valid1 = MAE(test_1[2],OS_time_valid1_knn)
r2_knn_valid1 = R2(test_1[2],OS_time_valid1_knn)
cor_knn_valid1 = cor(test_1[2],OS_time_valid1_knn)

test_results_valid1_knn<-cbind(cor_knn_valid1,rmse_knn_valid1,r2_knn_valid1,mae_knn_valid1)
write.csv(test_results_valid1_knn,file = "test_results_valid1_knn",row.names=F, quote = F)
#######################################################################################################


############################## Decision Tree (DT) Model ###################################
#################### Model Developement using svmRadial method at best parameters ################################
model_DT<- train(reg_function, data = train1 , method = "rpart", trControl = train.crossVal)

##################### save model #################################
save(model_DT, file = "vimp_DT_model.rda")

##################### extract best parameters of the final model #######################
parameters_DT<-model_DT$bestTune
write.csv(parameters_DT,file = "Best_DT_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
DT_rmse<-model_DT$results$RMSE
DT_r2<-model_DT$results$Rsquared
DT_mae<-model_DT$results$MAE
perf_DT<-cbind(DT_rmse,DT_r2,DT_mae)
write.csv(perf_DT,file = "perf_train_DT",row.names=F, quote = F)



##################### Prediction on Validation dataset1 ########################
test1_pred_DT<-predict(model_DT, newdata = valid_data1)
####################### round upto 2 digit ##############################
OS_time_valid1_DT<-round(test1_pred_DT,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_DT= cbind(valid_data1_id,OS_time_valid1_DT)
write.csv(Prediction_result_valid1_DT,file = "OS_time_Prediction_VH_vimp_valid1_DT",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_DT_valid1 = rmse(test_1[2],OS_time_valid1_DT)
mae_DT_valid1 = MAE(test_1[2],OS_time_valid1_DT)
r2_DT_valid1 = R2(test_1[2],OS_time_valid1_DT)
cor_DT_valid1 = cor(test_1[2],OS_time_valid1_DT)

test_results_valid1_DT<-cbind(cor_DT_valid1,rmse_DT_valid1,r2_DT_valid1,mae_DT_valid1)
#write.csv(test_results_valid1_DT,file = "test_results_valid1_DT",row.names=F, quote = F)
#######################################################################################################


############################## Ridge Model ###################################
#################### Model Developement using ridge method at best parameters ################################
lambda <- 10^seq(-3, 3, length = 100)
model_ridge<- train(reg_function, data = train1 , method = "glmnet", trControl = train.crossVal, tuneGrid = expand.grid(alpha = 0, lambda = lambda))

##################### save model #################################
save(model_ridge, file = "vimp_ridge_model.rda")

##################### extract best parameters of the final model #######################
parameters_ridge<-model_ridge$bestTune
write.csv(parameters_ridge,file = "Best_ridge_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
ridge_rmse<-model_ridge$results$RMSE
ridge_r2<-model_ridge$results$Rsquared
ridge_mae<-model_ridge$results$MAE
perf_ridge<-cbind(ridge_rmse,ridge_r2,ridge_mae)
write.csv(perf_ridge,file = "perf_train_ridge",row.names=F, quote = F)



##################### Prediction on Validation dataset1 ########################
test1_pred_ridge<-predict(model_ridge, newdata = valid_data1)
####################### round upto 2 digit ##############################
OS_time_valid1_ridge<-round(test1_pred_ridge,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_ridge= cbind(valid_data1_id,OS_time_valid1_ridge)
write.csv(Prediction_result_valid1_ridge,file = "OS_time_Prediction_VH_vimp_valid1_ridge",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_ridge_valid1 = rmse(test_1[2],OS_time_valid1_ridge)
mae_ridge_valid1 = MAE(test_1[2],OS_time_valid1_ridge)
r2_ridge_valid1 = R2(test_1[2],OS_time_valid1_ridge)
cor_ridge_valid1 = cor(test_1[2],OS_time_valid1_ridge)

test_results_valid1_ridge<-cbind(cor_ridge_valid1,rmse_ridge_valid1,r2_ridge_valid1,mae_ridge_valid1)
#write.csv(test_results_valid1_ridge,file = "test_results_valid1_ridge",row.names=F, quote = F)
#######################################################################################################



############################## Elastic Net (EN) Model ###################################
#################### Model Developement using EN method at best parameters ################################
lambda <- 10^seq(-3, 3, length = 100)
model_EN<- train(reg_function, data = train1 , method = "glmnet", trControl = train.crossVal, tuneGrid = expand.grid(alpha = 0.5, lambda = lambda))

##################### save model #################################
save(model_EN, file = "vimp_EN_model.rda")

##################### extract best parameters of the final model #######################
parameters_EN<-model_EN$bestTune
write.csv(parameters_EN,file = "Best_EN_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
EN_rmse<-model_EN$results$RMSE
EN_r2<-model_EN$results$Rsquared
EN_mae<-model_EN$results$MAE
perf_EN<-cbind(EN_rmse,EN_r2,EN_mae)
write.csv(perf_EN,file = "perf_train_EN",row.names=F, quote = F)


##################### Prediction on Validation dataset1 ########################
test1_pred_EN<-predict(model_EN, newdata = valid_data1)
####################### round upto 2 digit ##############################
OS_time_valid1_EN<-round(test1_pred_EN,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_EN= cbind(valid_data1_id,OS_time_valid1_EN)
write.csv(Prediction_result_valid1_EN,file = "OS_time_Prediction_VH_vimp_valid1_EN",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_EN_valid1 = rmse(test_1[2],OS_time_valid1_EN)
mae_EN_valid1 = MAE(test_1[2],OS_time_valid1_EN)
r2_EN_valid1 = R2(test_1[2],OS_time_valid1_EN)
cor_EN_valid1 = cor(test_1[2],OS_time_valid1_EN)

test_results_valid1_EN<-cbind(cor_EN_valid1,rmse_EN_valid1,r2_EN_valid1,mae_EN_valid1)
#write.csv(test_results_valid1_EN,file = "test_results_valid1_EN",row.names=F, quote = F)
#######################################################################################################




############################## Lasso Model ###################################
#################### Model Developement using Lasso method at best parameters ################################
lambda <- 10^seq(-3, 3, length = 100)
model_lasso<- train(reg_function, data = train1 , method = "glmnet", trControl = train.crossVal, tuneGrid = expand.grid(alpha = 1, lambda = lambda))

##################### save model #################################
save(model_lasso, file = "vimp_lasso_model.rda")

##################### extract best parameters of the final model #######################
parameters_lasso<-model_lasso$bestTune
write.csv(parameters_lasso,file = "Best_lasso_parameters",row.names=F, quote = F)
##################### extract performance on training set ##############################
lasso_rmse<-model_lasso$results$RMSE
lasso_r2<-model_lasso$results$Rsquared
lasso_mae<-model_lasso$results$MAE
perf_lasso<-cbind(lasso_rmse,lasso_r2,lasso_mae)
write.csv(perf_lasso,file = "perf_train_lasso",row.names=F, quote = F)



##################### Prediction on Validation dataset1 ########################
test1_pred_lasso<-predict(model_lasso, newdata = valid_data1)
####################### round upto 2 digit ##############################
OS_time_valid1_lasso<-round(test1_pred_lasso,2)
###################Combining Sample Ids with Prediction values ################### 
Prediction_result_valid1_lasso= cbind(valid_data1_id,OS_time_valid1_lasso)
write.csv(Prediction_result_valid1_lasso,file = "OS_time_Prediction_VH_vimp_valid1_lasso",row.names=F, quote = F)
######################## Calculation of Performance measures of the model ################################
rmse_lasso_valid1 = rmse(test_1[2],OS_time_valid1_lasso)
mae_lasso_valid1 = MAE(test_1[2],OS_time_valid1_lasso)
r2_lasso_valid1 = R2(test_1[2],OS_time_valid1_lasso)
cor_lasso_valid1 = cor(test_1[2],OS_time_valid1_lasso)

test_results_valid1_lasso<-cbind(cor_lasso_valid1,rmse_lasso_valid1,r2_lasso_valid1,mae_lasso_valid1)
write.csv(test_results_valid1_lasso,file = "test_results_valid1_lasso",row.names=F, quote = F)

#######################################################################################################



