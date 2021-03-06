#### AP vs PP samples

## https://dataaspirant.com/2017/01/19/support-vector-machine-classifier-implementation-r-caret-package/

get_best_result <- function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  return(best_result)
}

# Attach Packages
library(tidyverse)    # data manipulation and visualization
library(kernlab)      # SVM methodology
library(e1071)        # SVM methodology
library(ISLR)         # contains example data set "Khan"
library(RColorBrewer) # customized coloring of plots
library(caret)
library(glmnet)
### try random forest to classify AP and PP
df_sample_data_AP_sorted = readRDS("./df_sample_data_AP_sorted.Rds")
df_sample_data_PP_sorted = readRDS("./df_sample_data_PP_sorted.Rds")


df_sample_data_AP_PP_case_ctrl = rbind(df_sample_data_AP_sorted[,c(1,7,seq(11,312))],
                                        df_sample_data_PP_sorted[,c(1,7,seq(11,312))]) 

df_sample_data_AP_PP_case_ctrl = df_sample_data_AP_PP_case_ctrl[,seq(2,304)]

set.seed(5011)

cv_5 = trainControl(method = "cv", number = 10)

AP_PP_case_ctrl_elnet = train(
  class ~ ., data = df_sample_data_AP_PP_case_ctrl,
  method = "glmnet",
  trControl = cv_5,
  tuneLength = 10
)

br_el = get_best_result(AP_PP_case_ctrl_elnet)

X = model.matrix(class ~ ., df_sample_data_AP_PP_case_ctrl)[, -1]
y = df_sample_data_AP_PP_case_ctrl$class

fit_el_cv = cv.glmnet(X, y, alpha = 0.5, family = "binomial")

coeff1 = coef(fit_el_cv, s="lambda.min")
features1 = NULL
coeffs1 = NULL
features1 = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
coeffs1   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included

res1 = NULL
res1$features = features1[-1]
res1$coefs = coeffs1[-1]

id_list_1 = as.integer(unlist(lapply(res1$features, function(z){str_split(z,"_")[[1]][2]})))
df_compound = fread("./compound_list.csv")
AP_PP_compund_name_list = df_compound[id_list_1, 2]

saveRDS(id_list_1,"./id_list_AP_PP.Rds")


df_plot = df_sample_data_AP_PP[,id_list_1]
df_plot$type = c(rep('0',length(df_sample_data_AP_sorted$sample_id)),
                 rep('1',length(df_sample_data_PP_sorted$sample_id)))
df_plot$sample_id = c(df_sample_data_AP_sorted$sample_id, df_sample_data_PP_sorted$sample_id)

df_plot_melt = melt(df_plot, id=c("sample_id","type"))
names(df_plot_melt) = c("sample_id", "type", "compound", "value")   


ggplot(data=df_plot_melt, aes(x=compound, y=value, group=sample_id, color=type)) +
  geom_line() + geom_point() + ylab("log-value") +
  scale_color_manual(name="Type", labels=c("AP","PP"), values = c("#E69F00","#56B4E9")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()


######################
df_plot_1 = df_plot[-dim(df_plot)[2]]
set.seed(97333)
intrain = createDataPartition(y = df_plot$type, p = 0.8, list = FALSE)
training = df_plot_1[intrain,]
testing  = df_plot_1[-intrain,]

trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
set.seed(52133)
svm_Linear <- train(type ~., data = training, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 5)

test_pred <- predict(svm_Linear, newdata = testing)
confusionMatrix(factor(test_pred), factor(testing$type) )

grid <- expand.grid(C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5))
set.seed(3733)
svm_Linear_Grid <- train(type ~., data = training, method = "svmLinear",
                         trControl=trctrl,
                         preProcess = c("center", "scale"),
                         tuneGrid = grid,
                         tuneLength = 5)
plot(svm_Linear_Grid)
test_pred_grid <- predict(svm_Linear_Grid, newdata = testing)
confusionMatrix(factor(test_pred_grid), factor(testing$type) )

set.seed(3133)
svm_Radial <- train(type ~., data = training, method = "svmRadial",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 5)

test_pred_radial <- predict(svm_Radial, newdata = testing)
confusionMatrix(factor(test_pred_radial), factor(testing$type) )

grid_radial <- expand.grid(sigma = c(0.01, 0.02, 0.025, 0.03, 0.04,
                                     0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                           C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                                 1, 1.5, 2,5))
set.seed(3233)
svm_Radial_Grid <- train(type ~., data = training, method = "svmRadial",
                         trControl=trctrl,
                         preProcess = c("center", "scale"),
                         tuneGrid = grid_radial,
                         tuneLength = 5)
test_pred_grid_radial <- predict(svm_Radial_Grid, newdata = testing)
confusionMatrix(factor(test_pred_grid_radial), factor(testing$type) )

