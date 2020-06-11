##### stanford metabolomics data 

library(data.table)
library(shiny)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(glmnet)
library(gtools)

set.seed(96357)

setwd("~/Desktop/ToDo/CPG/GEMM_Alan/metabolomics/")

df_compound = fread("./compound_list.csv")
df_key = fread("./sample_key.csv")
df_sample = fread("./gemm_sample_info.txt")
df_subject = fread("./gemm_subject_info.txt")
#######################

N_DAYS_OF_WEEK = 7
### data preparation ###
col_id_start = 7
col_id_end   = 146

s_size = col_id_end - col_id_start + 1

L = dim(df_key)[1]

df_sample_data = NULL
### AP is 0 and PP is 1 and I consider days_wrt_delivery = 0 as PP
#### data not avilable for depress_46 and depress_88
for(i in 1:L){
  id = which(df_sample$RID == df_key$RID[i] & 
               df_sample$shipping_box == df_key$shipping_box[i] &
               df_sample$position_in_shipping_box == df_key$position_in_shipping_box[i]) 
  
  df_sample_data$sample_id[i] = as.character(df_key$sample_id[i])
  df_sample_data$RID[i] =  df_sample$RID[id]
  df_sample_data$days_wrt_delivery[i] = df_sample$days_wrt_delivery[id]
  df_sample_data$week_wrt_delivery[i] = df_sample$days_wrt_delivery[id]/N_DAYS_OF_WEEK
  df_sample_data$AP_or_PP[i] = as.integer(df_sample$days_wrt_delivery[id]>=0)
  df_sample_data$type[i] =  df_sample$type[id]
  df_sample_data$box[i] = df_sample$shipping_box[id]
  df_sample_data$pos[i] = df_sample$position_in_shipping_box[id]
  df_sample_data$data_available[i] = as.integer(!((df_sample_data$sample_id[i]=='depress_46') | 
                                                  (df_sample_data$sample_id[i]=='depress_88')))
  
}
df_sample_data =  as.data.frame(df_sample_data)

########### extract data from df_compound
#### total samples
X = t(df_compound[, col_id_start:col_id_end])
Y = unlist(lapply(as.vector(rownames(X)), function(z){id=which(df_sample_data$sample_id == z);
                                            return(df_sample_data$week_wrt_delivery[id])}))

### AP samples
AP_sample_id = as.vector(df_sample_data$sample_id[which(df_sample_data$AP_or_PP == 0)])
AP_row_id = unlist(lapply(AP_sample_id, function(z){ id = which(as.vector(rownames(X)) == z); return(id)}))

X1 = X[AP_row_id,]
Y1 = Y[AP_row_id]

### PP samples
PP_sample_id = as.vector(df_sample_data$sample_id[which(df_sample_data$AP_or_PP == 1)])
PP_row_id = unlist(lapply(PP_sample_id, function(z){ id = which(as.vector(rownames(X)) == z); return(id)}))

X2 = X[PP_row_id,]
Y2 = Y[PP_row_id]

###################s

####### glmnet on AP samples
  
  cat('dim X1: ', dim(X1)," length Y1: ", length(Y1),'\n')
  
  colnames(X1) = as.character(seq(1,dim(X1)[2]))
  
  AP_compund_name_list_arr = NULL
  
  a1 = 0.5
  i = 0
  for(a in seq(0,1.0,by=0.1)){
    i = i+1
    a1 = a
    
    fit1 = cv.glmnet(X1, Y1, nfolds = nrow(X1), alpha = a1, family = "gaussian", grouped = FALSE)
    
    coeff1 = coef(fit1, s="lambda.min")
    ###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
    features1 = NULL
    coeffs1 = NULL
    features1 = coeff1@Dimnames[[1]][ which(coeff1 != 0 ) ]  #intercept included
    coeffs1   = coeff1              [ which(coeff1 != 0 ) ]  #intercept included
    
    res1 = NULL
    res1$features = features1[-1]
    res1$coefs = coeffs1[-1]
    
    id_list_1 = as.integer(res1$features)
    AP_compund_name_list = df_compound[id_list_1, 2]
    
    AP_compund_name_list_arr[[i]] = AP_compund_name_list
    ####
}

###### glmnet on PP samples
    
    # cat('dim X2: ', dim(X2)," length Y2: ", length(Y2),'\n')
    # 
    # colnames(X2) = as.character(seq(1,dim(X2)[2]))    
    # a2 = 0.5
    # 
    # #for(a in seq(0,1.0,by=0.1)){
    # 
    # 
    # fit2 = cv.glmnet(X2, Y2, nfolds = nrow(X2), alpha = a2, family = "gaussian", grouped=FALSE)
    # 
    # coeff2 = coef(fit2, s="lambda.min")
    # ###https://stackoverflow.com/questions/27801130/extracting-coefficient-variable-names-from-glmnet-into-a-data-frame
    # features2 = NULL
    # coeffs2 = NULL
    # features2 = coeff2@Dimnames[[1]][ which(coeff2 != 0 ) ]  #intercept included
    # coeffs2   = coeff2              [ which(coeff2 != 0 ) ]  #intercept included
    # 
    # res2 = NULL
    # res2$features = features2[-1]
    # res2$coefs = coeffs2[-1]
    # 
    # id_list_2 = as.integer(res2$features)
    # PP_compund_name_list = df_compound[id_list_2, 2]
    # 
    # ############ 
    # 
    # 
    # 
