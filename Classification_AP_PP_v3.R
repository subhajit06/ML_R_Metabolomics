##### PLS_DA analysis

library(mixOmics)
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

df_compound = fread("./compound_list.csv")

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

df_sample_data_AP_PP = df_sample_data_AP_PP_1[, seq(3,305)]  ### 302 variates and y


X <- as.matrix(df_sample_data_AP_PP[,1:302])
Y <- ifelse((df_sample_data_AP_PP$y == 'yes') , 1, 0)
Y = as.factor(df_sample_data_AP_PP$y)

ppmd.plsda <- plsda(X, Y, ncomp = 20)  # set ncomp to 10 for performance assessment later
plotIndiv(ppmd.plsda , comp = 1:2,
          group = df_sample_data_AP_PP$y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'PLSDA on SRBCT')

list.keepX <- c(1:10,  seq(10, 50, 2))
tune.splsda.ppmd <- tune.splsda(X, Y, ncomp = 6, validation = 'Mfold', folds = 5, 
                                 progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                 test.keepX = list.keepX, nrepeat = 10, cpus = 2)

splsda.ppmd = splsda(X, Y, ncomp = tune.splsda.ppmd$choice.ncomp$ncomp, 
                     keepX = tune.splsda.ppmd$choice.keepX[1:ncomp]) 

plotIndiv(splsda.ppmd, comp = c(1,2),
          group = Y, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'sPLS-DA on SRBCT, comp 1 & 2')

auroc(splsda.ppmd, roc.comp = tune.splsda.ppmd$choice.ncomp$ncomp)
set.seed(40) # for reproducibility, only when the `cpus' argument is not used
# takes about 1 min to run
perf.ppmd <- perf(splsda.ppmd, validation = "Mfold", folds = 5,
                   dist = 'max.dist', nrepeat = 10,
                   progressBar = FALSE) 
perf.ppmd$error.rate
plot(perf.ppmd, col = color.mixo(5))


# here we match the selected variables to the stable features
ind.match = match(selectVar(splsda.ppmd, comp = 1)$name, 
                  names(perf.ppmd$features$stable[[1]]))
#extract the frequency of selection of those selected variables
Freq = as.numeric(perf.ppmd$features$stable[[1]][ind.match])

data.frame(selectVar(splsda.ppmd, comp = 1)$value, Freq)

id_list_AP_PP_new = readRDS("./id_list_AP_PP_new.Rds")
AP_PP_compund_name_list = df_compound[id_list_AP_PP_new, 2]

id1 = as.integer(unlist(lapply(selectVar(splsda.ppmd, comp = 1)$name,function(x){
  str_split(x,"comp_")[[1]][2]})))

comp_1_compund_list = df_compound$Compound.name[id1]
intersect(id_list_AP_PP_new,id1)

id2 = as.integer(unlist(lapply(selectVar(splsda.ppmd, comp = 2)$name,function(x){
  str_split(x,"comp_")[[1]][2]})))

comp_2_compund_list = df_compound$Compound.name[id2]

intersect(id_list_AP_PP_new,id2)

pdf("./plots/cim_ppmd_comp_1.pdf")
cim(splsda.ppmd, comp=1, title ="Component 1")
dev.off()

pdf("./plots/arrow_ppmd.pdf")
plotArrow(splsda.ppmd, legend=T)
dev.off()
length(intersect(id_list_AP_PP_new,union(id1,id2)))

plotLoadings(splsda.ppmd, comp = 1, title = 'Loadings on comp 1', 
             contrib = 'max', method = 'mean')

plotLoadings(splsda.ppmd, comp = 2, title = 'Loadings on comp 2', 
             contrib = 'max', method = 'mean')

