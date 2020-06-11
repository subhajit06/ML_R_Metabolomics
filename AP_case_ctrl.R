#### AP: case vs control samples

## https://daviddalpiaz.github.io/r4sl/elastic-net.html#regression-1
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

df_sample_data_AP_sorted = readRDS("./df_sample_data_AP_sorted.Rds")

df_sample_data_AP = df_sample_data_AP_sorted[, c(7,seq(11,312))]

set.seed(92713)

cv_5 = trainControl(method = "cv", number = 5)

AP_case_ctrl_elnet = train(
  class ~ ., data = df_sample_data_AP,
  method = "glmnet",
  trControl = cv_5,
  tuneLength = 5
)

br_el = get_best_result(AP_case_ctrl_elnet)

X = model.matrix(class ~ ., df_sample_data_AP)[, -1]
y = df_sample_data_AP$class

fit_el_cv = cv.glmnet(X, y, alpha = br_el$alpha, family = "binomial")

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
AP_case_ctrl_compund_name_list = df_compound[id_list_1, 2]

saveRDS(id_list_1,"./id_list_AP_case_ctrl.Rds")


df_plot = df_sample_data_AP[,id_list_1]
df_plot$type = df_sample_data_AP$class
df_plot$sample_id = df_sample_data_AP_sorted$sample_id

df_plot_melt = melt(df_plot, id=c("sample_id","type"))
names(df_plot_melt) = c("sample_id", "type", "compound", "value")   

ggplot(data=df_plot_melt, aes(x=compound, y=value, group=sample_id, color=type)) +
  geom_line() + geom_point() + ylab("log-value") + 
  scale_color_manual(name="Type", labels=c("Case","Control"), values = c("#E69F00","#56B4E9")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()


######################
df_plot_1 = df_plot[-dim(df_plot)[2]]
id1 = which(df_plot_1$type == 'No')
df_plot_1$type[id1] = 0
df_plot_1$type[-id1] = 1
df_plot_1$type = as.factor(df_plot_1$type)
  
set.seed(65337)
intrain = createDataPartition(y = df_plot_1$type, p = 0.75, list = FALSE)
training = df_plot_1[intrain,]
testing  = df_plot_1[-intrain,]

trctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
set.seed(89133)
svm_Linear <- train(type ~., data = training, method = "svmLinear",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 5)

test_pred <- predict(svm_Linear, newdata = testing)
confusionMatrix(factor(test_pred), factor(testing$type) )

grid <- expand.grid(C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5))
svm_Linear_Grid <- train(type ~., data = training, method = "svmLinear",
                         trControl=trctrl,
                         preProcess = c("center", "scale"),
                         tuneGrid = grid,
                         tuneLength = 5)
plot(svm_Linear_Grid)
test_pred_grid <- predict(svm_Linear_Grid, newdata = testing)
confusionMatrix(factor(test_pred_grid), factor(testing$type) )

svm_Radial <- train(type ~., data = training, method = "svmRadial",
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 5)

test_pred_radial <- predict(svm_Radial, newdata = testing)
confusionMatrix(factor(test_pred_radial), testing$type )


