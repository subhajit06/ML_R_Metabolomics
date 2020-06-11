#### AP vs PP samples
#### removing Delivery samples
### http://chenyuan.date/2017/11/Penalized-Regression-Ridge-Lasso-Elastic-Net/
### https://daviddalpiaz.github.io/r4sl/elastic-net.html

#rm(list=ls()) 

# Attach Packages
library(tidyverse)    # data manipulation and visualization
library(kernlab)      # SVM methodology
library(e1071)        # SVM methodology
library(data.table)
library(RColorBrewer) # customized coloring of plots
library(caret)
library(glmnet)
library(gridExtra)
library(grid)
library(corrplot)
library(factoextra)
library(ggfortify)

calc_accuracy = function(actual, predicted) {
  mean(actual == predicted)
}

get_best_result <- function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  return(best_result)
}

### try random forest to classify AP and PP (DL are technically PP data)
df_sample_data_AP_sorted = readRDS("./df_sample_data_AP_sorted.Rds")
df_sample_data_PP_sorted = readRDS("./df_sample_data_PP_sorted.Rds")
#df_sample_data_DL_sorted = readRDS("./df_sample_data_DL_sorted.Rds")

#df_sample_data_AP_PP_1 = rbind(df_sample_data_AP_sorted[,c(1,5,seq(11,312))],
#                                df_sample_data_PP_sorted[,c(1,5,seq(11,312))])

#### include 0 days as PP samples
df_sample_data_AP_PP_1 = rbind(df_sample_data_AP_sorted[,c(1,5,seq(11,312))],
                               #df_sample_data_DL_sorted[,c(1,5,seq(11,312))],
                               df_sample_data_PP_sorted[,c(1,5,seq(11,312))]) 

df_sample_data_AP_PP_1$y = rep('yes', dim(df_sample_data_AP_PP_1)[1])
df_sample_data_AP_PP_1$y[1:dim(df_sample_data_AP_sorted)[1]] = rep('no', dim(df_sample_data_AP_sorted)[1])

df_sample_data_AP_PP = df_sample_data_AP_PP_1[, seq(3, 305)]  ### 302 variates and y

#set.seed(87317)
#set.seed(87319)
#set.seed(371319)
#set.seed(24387)
#set.seed(75513)
set.seed(873175)

########### training and testing data split by y which is the response variable
training <- df_sample_data_AP_PP$y %>%  createDataPartition(p = 0.8, list = FALSE)

train_data <- df_sample_data_AP_PP[training, ]
test_data  <- df_sample_data_AP_PP[-training, ]

## Predictor variables
#X_AP_PP <- model.matrix(y~., train_data)[,-1]  
## Response variable
#y_AP_PP <- ifelse(train_data$y == "yes", 1, 0)


cv_4 = trainControl(method = "cv", number = 4)
cv_r_4 = trainControl(method = "repeatedcv", number = 4, repeats = 3)

en_AP_PP <- train(
  y ~., data = train_data, 
  method = "glmnet",
  family="binomial",
  #relax = TRUE,
  trControl = cv_r_4,
  tuneLength = 4
)

br_el = get_best_result(en_AP_PP)

# Best tuning parameter
en_AP_PP$bestTune
coef(en_AP_PP$finalModel, en_AP_PP$bestTune$lambda)

################ testing
#x.test <- model.matrix(y ~., test_data)[,-1]
#predictions <- en_AP_PP %>% predict(x.test)

predicted = predict(en_AP_PP, newdata = test_data)
acc_test = calc_accuracy(actual = test_data$y, predicted)
cm = confusionMatrix(predicted, as.factor(test_data$y))

#which(coef(en_AP_PP$finalModel, en_AP_PP$bestTune$lambda) != 0)

### 153 162 205 274 281

########

coeff1 = coef(en_AP_PP$finalModel, en_AP_PP$bestTune$lambda)
features1 = NULL
coeffs1 = NULL
features1 = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
coeffs1   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included

res1 = NULL
res1$features = features1[-1]
res1$coefs = coeffs1[-1]

id_list_AP_PP = as.integer(unlist(lapply(res1$features, function(z){str_split(z,"_")[[1]][2]})))
df_compound = fread("./compound_list.csv")
AP_PP_compund_name_list = df_compound[id_list_AP_PP, 2]
saveRDS(id_list_AP_PP,"./id_list_AP_PP_new_1.Rds")


id_list_AP_PP_old = readRDS("./id_list_AP_PP.Rds")
id_list_AP_PP_new_0 = readRDS("./id_list_AP_PP_new.Rds")

AP_PP_compund_name_list_old = df_compound[id_list_AP_PP_old, 2]
AP_PP_compund_name_list_new_0 = df_compound[id_list_AP_PP_new_0, 2]

intersect(AP_PP_compund_name_list_old, AP_PP_compund_name_list)
intersect(AP_PP_compund_name_list_new_0, AP_PP_compund_name_list)

id_list_AP_PP
AP_PP_compund_name_list

######################
######################
pdf("./plots/compound_table_AP_PP_classify_new_1.pdf", width=12, height=5)
grid.table(df_compound[id_list_AP_PP,1:6])
dev.off()

id_list_intersect = intersect(id_list_AP_PP_new_0,id_list_AP_PP)
pdf("./plots/compound_table_AP_PP_classify_intersect.pdf", width=12, height=5)
grid.table(df_compound[id_list_intersect,1:6])
dev.off()


df_plot = df_sample_data_AP_PP[,id_list_AP_PP]
df_plot$type = c(rep('0',length(df_sample_data_AP_sorted$sample_id)),
                 #rep('1',length(df_sample_data_DL_sorted$sample_id)),
                 rep('1',length(df_sample_data_PP_sorted$sample_id)))

df_plot$sample_id = c(df_sample_data_AP_sorted$sample_id, 
                      #df_sample_data_DL_sorted$sample_id,
                      df_sample_data_PP_sorted$sample_id)
df_plot = as.data.frame(df_plot)

df_gemm_selected_features = df_plot
saveRDS(df_gemm_selected_features, "./df_gemm_selected_features_1.Rds")

df_plot_melt = reshape2::melt(df_plot, id=c("sample_id","type"))
names(df_plot_melt) = c("sample_id", "type", "compound", "value")   

pdf("./plots/AP_PP_sample_classify_from_metabolites_new_1.pdf",width=12, height=8)
ggplot(data=df_plot_melt, aes(x=compound, y=value, group=sample_id, color=type)) +
  geom_line() + geom_point() + ylab("log-value") +
  scale_color_manual(name="Type", labels=c("AP","PP"), values = c("#E69F00","#56B4E9")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()
dev.off()  

################### create a matrix with all selected features for 134 samples 
comp_id_arr = unlist(lapply(res1$features, function(x) {which(colnames(df_sample_data_AP_PP_1) == x)}))
df_sample_data_selected_comp = df_sample_data_AP_PP_1[, c(1,2,comp_id_arr)]
saveRDS(df_sample_data_selected_comp,"./df_sample_data_selected_comp_1.Rds")

######## PCA plot with those 13 features
df_AP_PP.pr <- prcomp(df_sample_data_selected_comp[c(3:13)], center = TRUE, scale = TRUE)
summary(df_AP_PP.pr)

#autoplot(df_AP_PP.pr, data = df_sample_data_selected_comp, colour = 'AP_or_PP')

pdf("./plots/scree_plot_13_features_1.pdf",width = 10, height=8)
fviz_eig(df_AP_PP.pr)
dev.off()

groups = ifelse(df_sample_data_selected_comp$AP_or_PP == 0, 'AP','PP')
pdf("./plots/2D_pca_plot_13_features_1.pdf", width = 10, height=8)
fviz_pca_ind(df_AP_PP.pr,
             col.ind = groups, # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             palette = c("blue", "red"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             geom = "point",
             pointshape = 20, 
             pointsize = 2,
             repel = TRUE)+
  ggtitle("2D PCA-plot from 11 feature dataset") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


###############

library(data.table)
library(ggplot2)
library(tidyverse)
library(caret)
library(dplyr)         # Used by caret
library(kernlab)       # support vector machine 
library(pROC)	       # plot the ROC curves
library(randomForest)

setwd("~/Desktop/ToDo/CPG/GEMM_Alan/metabolomics/")

df_compound = fread("./compound_list.csv")
id_list_AP_PP_1 = readRDS("./id_list_AP_PP_new_1.Rds")
df_gemm_selected_features = readRDS("./df_gemm_selected_features_1.Rds")

df_gemm = df_gemm_selected_features %>% dplyr::select(!"sample_id")

df_gemm$type[df_gemm$type==1] = "yes"
df_gemm$type[df_gemm$type==0] = "no"
df_gemm$type = as.factor(df_gemm$type)

set.seed(15392)

#df_compound[id_list_AP_PP,1:6]

### Get the Data
# Load the data and construct indices to divide it into training and test data sets.
trainIndex <- createDataPartition(df_gemm$type, p=0.8,list=FALSE)
trainData <- df_gemm[trainIndex,]
testData  <- df_gemm[-trainIndex,]
trainX <-trainData[,1:11]        # Pull out the variables for training
sapply(trainX,summary)           # Look at a summary of the training data

## SUPPORT VECTOR MACHINE MODEL


################# Linear SVM model
# First pass

# Setup for cross validation
ctrl <- trainControl(method="repeatedcv",   # 10fold cross validation
                     repeats=5,		    # do 5 repititions of cv
                     #summaryFunction=twoClassSummary,	# Use AUC to pick the best model
                     classProbs=TRUE)
set.seed(15392)
svmLinear <- train(type ~ .,
                   data = trainData, 
                   method = "svmLinear",
                   preProcess = c("center","scale"),
                   trControl = ctrl,
                   tuneLength = 5)


set.seed(15397)
svmLinear.grid <- train(type ~ .,
                        data = trainData, 
                        method = "svmLinear",
                        trControl = ctrl,
                        tuneLength = 5,
                        tuneGrid = expand.grid(C = seq(0.1, 2, length = 10)),
                        preProcess = c("center","scale"))
svmLinear
svmLinear.grid

predicted.classes.svmLinear <- svmLinear %>% predict(testData)
accuracy_svmLinear = mean(predicted.classes.svmLinear == testData$type)
cmSVMLinear <-confusionMatrix(predicted.classes.svmLinear, testData$type)
importance.svmLinear <- varImp(svmLinear, scale=FALSE)
plot(importance.svmLinear)

predicted.classes.svmLinearGrid <- svmLinear.grid %>% predict(testData)
accuracy_svmLinearGrid = mean(predicted.classes.svmLinearGrid == testData$type)
cmSVMLinearGrid <-confusionMatrix(predicted.classes.svmLinearGrid, testData$type)
importance.svmLinearGrid <- varImp(svmLinear.grid, scale=FALSE)
plot(importance.svmLinearGrid)

#### random forest
set.seed(13782)

ctrl = trainControl(method="repeatedcv",   # 10fold cross validation
                    repeats=5)

rf_model = train(type ~ .,
                 data=trainData,
                 method="rf",
                 trControl=ctrl,
                 tuneLength = 5,
                 #tuneGrid=rf_grid,
                 prox=TRUE)

predicted.classes.rf = rf_model %>% predict(testData) 

accuracy_rf = mean(predicted.classes.rf == testData$type)
cmRF <- confusionMatrix(predicted.classes.rf, testData$type)

set.seed(13782)
#mtry <- sqrt(ncol(trainX))
#tunegrid <- expand.grid(.mtry=mtry)
rf_grid <- expand.grid(mtry = c(1, 2, 3, 4, 5, 6))

rf_model.grid = train(type ~ .,
                      data=trainData,
                      method="rf",
                      trControl=ctrl,
                      tuneLength = 5,
                      tuneGrid=rf_grid
                      #prox=TRUE
)

predicted.classes.rf.grid = rf_model.grid %>% predict(testData) 

accuracy_rf.grid = mean(predicted.classes.rf.grid == testData$type)
cmRFGrid <- confusionMatrix(predicted.classes.rf.grid, testData$type)


