############ data description plot
library(data.table)
library(ggplot2)
library(dplyr)
library(caret)
library(glmnet)
library(gridExtra)
library(grid)

get_best_result <- function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

df_compound = fread("./compound_list.csv")
df_sample_data_AP_sorted = readRDS("./df_sample_data_AP_sorted.Rds")
df_sample_data_PP_sorted = readRDS("./df_sample_data_PP_sorted.Rds")

df_sample_data_AP_sorted_filtered = df_sample_data_AP_sorted[,c(4,seq(11,312))]
df_sample_data_PP_sorted_filtered = df_sample_data_PP_sorted[,c(4,seq(11,312))]

############### AP samples

set.seed(13057)

cv_4 = trainControl(method = "cv", number = 4)

AP_GA_elnet = train(
  week_wrt_delivery ~ ., data = df_sample_data_AP_sorted_filtered,
  method = "glmnet",
  trControl = cv_4,
  tuneLength = 4
)

br_el = get_best_result(AP_GA_elnet)

X = model.matrix(week_wrt_delivery ~ ., df_sample_data_AP_sorted_filtered)[, -1]
y = df_sample_data_AP_sorted_filtered$week_wrt_delivery

fit_el_cv = cv.glmnet(X, y, alpha = br_el$alpha, family = "gaussian", type.measure="mse")
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
saveRDS(AP_compund_name_list_1, './AP_compund_name_list_3_weeks.Rds')
AP_compund_name_list_1


pdf("./plots/compound_table_AP_GA.pdf", width=11.5, height=8)
grid.table(df_compound[id_list_AP,1:6])
dev.off()

id_list = unlist(lapply(id_list_AP, function(x){ 
  which(colnames(df_sample_data_AP_sorted) == paste0('comp_',x))}))

df_plot = df_sample_data_AP_sorted[,c(1,4,id_list)]

df_plot_melt = melt(df_plot, id=c("sample_id","week_wrt_delivery"))
names(df_plot_melt) = c("sample_id", "week_wrt_delivery", "compound", "value")   

pdf("./plots/AP_sample_GA_regression.pdf",width=12, height=8)
ggplot(data=df_plot_melt, aes(x=week_wrt_delivery, y=value, color = compound)) +
  geom_line() + geom_point() + ylab("log-value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()
dev.off()

############## PP samples

set.seed(74579)

cv_4 = trainControl(method = "cv", number = 4)
PP_GA_elnet = train(
  week_wrt_delivery ~ ., data = df_sample_data_PP_sorted_filtered,
  method = "glmnet",
  trControl = cv_4,
  tuneLength = 4
)

br_el = get_best_result(PP_GA_elnet)

X = model.matrix(week_wrt_delivery ~ ., df_sample_data_PP_sorted_filtered)[, -1]
y = df_sample_data_PP_sorted_filtered$week_wrt_delivery

fit_el_cv = cv.glmnet(X, y, alpha = br_el$alpha, family = "gaussian", type.measure="mse")
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
saveRDS(PP_compund_name_list_1, './PP_compund_name_list_3_weeks.Rds')
PP_compund_name_list_1

pdf("./plots/compound_table_PP_GA.pdf", width=11.5, height=10)
grid.table(df_compound[id_list_PP,1:6])
dev.off()


id_list = unlist(lapply(id_list_PP, function(x){ 
  which(colnames(df_sample_data_PP_sorted) == paste0('comp_',x))}))

df_plot = df_sample_data_PP_sorted[,c(1,4,id_list)]

df_plot_melt = melt(df_plot, id=c("sample_id","week_wrt_delivery"))
names(df_plot_melt) = c("sample_id", "week_wrt_delivery", "compound", "value")   

pdf("./plots/PP_sample_GA_regression.pdf",width=12, height=8)
ggplot(data=df_plot_melt, aes(x=week_wrt_delivery, y=value, color=compound)) +
  geom_line() + geom_point() + ylab("log-value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()
dev.off()


