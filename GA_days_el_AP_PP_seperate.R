############ data description plot
library(data.table)
library(ggplot2)
library(dplyr)
library(caret)
library(glmnet)

get_best_result <- function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

df_compound = fread("./compound_list.csv")
df_sample_data_AP_sorted = readRDS("./df_sample_data_AP_sorted.Rds")
df_sample_data_PP_sorted = readRDS("./df_sample_data_PP_sorted.Rds")

df_sample_data_AP_sorted_filtered = df_sample_data_AP_sorted[,c(3,seq(11,312))]
df_sample_data_PP_sorted_filtered = df_sample_data_PP_sorted[,c(3,seq(11,312))]

###############

set.seed(139057)

cv_4 = trainControl(method = "cv", number = 4)

df_sample_data_AP_sorted_filtered$days_wrt_delivery = 
    -(df_sample_data_AP_sorted_filtered$days_wrt_delivery)
  
AP_GA_elnet = train(
  days_wrt_delivery ~ ., data = df_sample_data_AP_sorted_filtered,
  method = "glmnet",
  family = "poisson",
  trControl = cv_4,
  tuneLength = 4
)

br_el = get_best_result(AP_GA_elnet)

X = model.matrix(days_wrt_delivery ~ ., df_sample_data_AP_sorted_filtered)[, -1]
y = df_sample_data_AP_sorted_filtered$days_wrt_delivery

fit_el_cv = cv.glmnet(X, y, alpha = br_el$alpha, family = "poisson",type.measure="mse")
#sqrt(fit_el_cv$cvm[fit_el_cv$lambda == fit_el_cv$lambda.min]) # CV-RMSE minimum

coeff1 = coef(fit_el_cv, s="lambda.min")
###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
features1 = NULL
coeffs1 = NULL
features1 = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
coeffs1   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included

res1 = NULL
res1$features = features1[-1]
res1$coefs = coeffs1[-1]

id_list_AP = as.integer(unlist(lapply(res1$features, function(z){str_split(z,"_")[[1]][2]})))
AP_compund_name_list_1 = df_compound[id_list_AP, 2]
saveRDS(AP_compund_name_list_1, './AP_compund_name_list_poisson_days.Rds')
AP_compund_name_list_1

id_list = unlist(lapply(id_list_AP, function(x){ 
  which(colnames(df_sample_data_AP_sorted) == paste0('comp_',x))}))

df_plot = df_sample_data_AP_sorted[,c(1,3,id_list)]

df_plot_melt = melt(df_plot, id=c("sample_id","days_wrt_delivery"))
names(df_plot_melt) = c("sample_id", "days_wrt_delivery", "compound", "value")   

ggplot(data=df_plot_melt, aes(x=days_wrt_delivery, y=value, color = compound)) +
  geom_line() + geom_point() + #ylab("log-value") +
  #scale_color_manual(name="Type", labels=c("AP","PP"), values = c("#E69F00","#56B4E9")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()


##############

set.seed(74579)

cv_4 = trainControl(method = "cv", number = 4)
PP_GA_elnet = train(
  days_wrt_delivery ~ ., data = df_sample_data_PP_sorted_filtered,
  method = "glmnet",
  family = "poisson",
  trControl = cv_4,
  tuneLength = 4
)

br_el = get_best_result(PP_GA_elnet)

X = model.matrix(days_wrt_delivery ~ ., df_sample_data_PP_sorted_filtered)[, -1]
y = df_sample_data_PP_sorted_filtered$days_wrt_delivery

fit_el_cv = cv.glmnet(X, y, alpha = br_el$alpha, family = "poisson", type.measure="mse")
#sqrt(fit_el_cv$cvm[fit_el_cv$lambda == fit_el_cv$lambda.min]) # CV-RMSE minimum

coeff1 = coef(fit_el_cv, s="lambda.min")
###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
features1 = NULL
coeffs1 = NULL
features1 = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
coeffs1   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included

res1 = NULL
res1$features = features1[-1]
res1$coefs = coeffs1[-1]

id_list_PP = as.integer(unlist(lapply(res1$features, function(z){str_split(z,"_")[[1]][2]})))
PP_compund_name_list_1 = df_compound[id_list_PP, 2]
saveRDS(PP_compund_name_list_1, './PP_compund_name_list_poisson_days.Rds')
PP_compund_name_list_1

id_list = unlist(lapply(id_list_PP, function(x){ 
  which(colnames(df_sample_data_PP_sorted) == paste0('comp_',x))}))

df_plot = df_sample_data_PP_sorted[,c(1,3,id_list)]

df_plot_melt = melt(df_plot, id=c("sample_id","days_wrt_delivery"))
names(df_plot_melt) = c("sample_id", "days_wrt_delivery", "compound", "value")   

ggplot(data=df_plot_melt, aes(x=days_wrt_delivery, y=value, color=compound)) +
  geom_line() + geom_point() + #ylab("log-value") +
  #scale_color_manual(name="Type", labels=c("AP","PP"), values = c("#E69F00","#56B4E9")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()



