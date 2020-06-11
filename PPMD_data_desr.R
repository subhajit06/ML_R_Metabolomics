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

df_sample_data = fread("./df_sample_data.tsv")
L = dim(df_sample_data)[2]

### AP samples
AP_row_id = which(df_sample_data$AP_or_PP == 0 & df_sample_data$days_wrt_delivery != 0 
                  & df_sample_data$data_available == 1)

df_compound = fread("./compound_list.csv")
col_id_start = 7
col_id_end   = 146
X = t(df_compound[, col_id_start:col_id_end])
colnames(X) = unlist(lapply(seq(1,dim(X)[2]), function(z){paste0('comp_',z)} ))

df_sample_data_AP = as.data.frame(matrix(NA, ncol = dim(X)[2]+dim(df_sample_data)[2], 
                                         nrow = length(AP_row_id)))

df_sample_data_AP[, 1:L] = df_sample_data[AP_row_id,]

row_select_AP = unlist(lapply(df_sample_data_AP$V1, function(z){which(rownames(X) == z)}))
X_AP = X[row_select_AP,]

df_sample_data_AP[,(dim(df_sample_data)[2]+1):dim(df_sample_data_AP)[2]] = X_AP

colnames(df_sample_data_AP)[1:L] = colnames(df_sample_data)
colnames(df_sample_data_AP)[(L+1):dim(df_sample_data_AP)[2]] = colnames(X_AP)

tmp_id = order(df_sample_data_AP$week_wrt_delivery)
df_sample_data_AP_sorted = df_sample_data_AP[tmp_id,]
saveRDS(df_sample_data_AP_sorted, "./df_sample_data_AP_sorted.Rds")

######## PP samples
PP_row_id = which(df_sample_data$AP_or_PP == 1 & df_sample_data$days_wrt_delivery != 0 
                  & df_sample_data$data_available == 1)

df_sample_data_PP = as.data.frame(matrix(NA, ncol = dim(X)[2]+dim(df_sample_data)[2], 
                                         nrow = length(PP_row_id)))

df_sample_data_PP[,1:L] = df_sample_data[PP_row_id,]

row_select_PP = unlist(lapply(df_sample_data_PP$V1, function(z){which(rownames(X) == z)}))
X_PP = X[row_select_PP,]

df_sample_data_PP[, (dim(df_sample_data)[2]+1):dim(df_sample_data_PP)[2]] = X_PP

colnames(df_sample_data_PP)[1:L] = colnames(df_sample_data)
colnames(df_sample_data_PP)[(L+1):dim(df_sample_data_PP)[2]] = colnames(X_PP)

tmp_id = order(df_sample_data_PP$week_wrt_delivery)
df_sample_data_PP_sorted = df_sample_data_PP[tmp_id,]

saveRDS(df_sample_data_PP_sorted, "./df_sample_data_PP_sorted.Rds")

############### DL samples

DL_row_id = which(df_sample_data$AP_or_PP == 1 & df_sample_data$days_wrt_delivery == 0 
                  & df_sample_data$data_available == 1)

df_compound = fread("./compound_list.csv")
col_id_start = 7
col_id_end   = 146
X = t(df_compound[, col_id_start:col_id_end])
colnames(X) = unlist(lapply(seq(1,dim(X)[2]), function(z){paste0('comp_',z)} ))

df_sample_data_DL = as.data.frame(matrix(NA, ncol = dim(X)[2]+dim(df_sample_data)[2], 
                                         nrow = length(DL_row_id)))

df_sample_data_DL[, 1:L] = df_sample_data[DL_row_id,]

row_select_DL = unlist(lapply(df_sample_data_DL$V1, function(z){which(rownames(X) == z)}))
X_DL = X[row_select_DL,]

df_sample_data_DL[,(dim(df_sample_data)[2]+1):dim(df_sample_data_DL)[2]] = X_DL

colnames(df_sample_data_DL)[1:L] = colnames(df_sample_data)
colnames(df_sample_data_DL)[(L+1):dim(df_sample_data_DL)[2]] = colnames(X_DL)

tmp_id = order(df_sample_data_DL$week_wrt_delivery)
df_sample_data_DL_sorted = df_sample_data_DL[tmp_id,]
saveRDS(df_sample_data_DL_sorted, "./df_sample_data_DL_sorted.Rds")



###############
### CV for AP data calculation

df_sample_data_AP_sorted_select = df_sample_data_AP_sorted[,c(4,seq(11,312))]

df_sample_data_AP_sorted_select_with_cv = 
  df_sample_data_AP_sorted_select %>% group_by(week_wrt_delivery) %>% 
  summarise_all(function(x) {mean(x)})

###############

df_X1_Y1 = as.data.frame(df_sample_data_AP_sorted_select_with_cv)

set.seed(3057)

cv_5 = trainControl(method = "cv", number = 5)

AP_GA_elnet = train(
  week_wrt_delivery ~ ., data = df_X1_Y1,
  method = "glmnet",
  trControl = cv_5,
  tuneLength = 5
)

br_el = get_best_result(AP_GA_elnet)

X = model.matrix(week_wrt_delivery ~ ., df_X1_Y1)[, -1]
y = df_X1_Y1$week_wrt_delivery

fit_el_cv = cv.glmnet(X, y, alpha = br_el$alpha, family = "gaussian")
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

id_list_1 = as.integer(unlist(lapply(res1$features, function(z){str_split(z,"_")[[1]][2]})))
AP_compund_name_list_1 = df_compound[id_list_1, 2]
saveRDS(AP_compund_name_list_1, './AP_compund_name_list_1.Rds')


##############

ggplot(df_sample_data_AP_sorted, aes(x = week_wrt_delivery, y = comp_1, group=1)) + 
  geom_line() + geom_point()





