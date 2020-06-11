library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(stringr)

df_sample_data = fread("./df_sample_data.tsv")
U = unique(df_sample_data$RID)
map_tbl = NULL
for(i in 1:length(U)){
  map_tbl$rid[i] = U[i]
  map_tbl$s_id[i] = sprintf('s_%d',i)
}
map_tbl = as.data.frame(map_tbl)

df_sample_data$s_id = unlist(lapply(df_sample_data$RID, function(x){
                                         id = which(map_tbl$rid == x); return(map_tbl$s_id[id])}))


tbl_rid = as.data.frame(table(df_sample_data$s_id))
tbl_rid_AP_PP = as.data.frame(with(df_sample_data, table(s_id, AP_or_PP)))

id1 = which(tbl_rid_AP_PP$AP_or_PP == 0)
id2 = which(tbl_rid_AP_PP$AP_or_PP == 1)

tbl_rid_AP_PP$timestamp[id1] = rep('AP', length(id1))
tbl_rid_AP_PP$timestamp[id2] = rep('PP', length(id2))

pdf("./plots/freq_samples.pdf")
ggplot(tbl_rid, aes(x=Freq)) + 
  geom_bar(stat = "count") + 
  scale_x_discrete(name ="Sample count", limits=c("1","2","3","4","5","6","7","8")) + 
  scale_y_discrete(name ="Number of subjects", limits=seq(1,20))
dev.off()

pdf("./plots/freq_samples_AP_PP.pdf",width=11,height=6)
ggplot(tbl_rid_AP_PP, aes(s_id, Freq, fill = timestamp)) + 
  geom_bar(stat="identity", width=0.7,position=position_dodge(0.5)) + 
  scale_fill_brewer(palette = "Dark2") + xlab("GEMM Sample ID") + ylab("Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

tbl_rid_AP_PP$Case = unlist(lapply(tbl_rid_AP_PP$s_id, function(x){
  df_sample_data$type[which(df_sample_data$s_id == x)[1]]}))

id1 = which(tbl_rid_AP_PP$Case == 'case')
id2 = which(tbl_rid_AP_PP$Case == 'control')

tbl_rid_AP_PP$Case[id1] = 'red'
tbl_rid_AP_PP$Case[id2] = 'blue'

pdf("./plots/freq_samples_AP_PP_case_ctrl.pdf",width=11,height=6)
ggplot(tbl_rid_AP_PP, aes(s_id, Freq, fill = timestamp)) + 
  geom_bar(stat="identity", width=0.7,position=position_dodge(0.5)) + 
  scale_fill_brewer(palette = "Dark2") + xlab("GEMM Sample ID") + ylab("Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, colour=tbl_rid_AP_PP$Case))
dev.off()


#### tabular view
length(which(df_sample_data$AP_or_PP == 0 & df_sample_data$data_available ==1))
length(which(df_sample_data$AP_or_PP == 1 & df_sample_data$data_available ==1))

length(which(df_sample_data$AP_or_PP == 0 & df_sample_data$data_available ==1 & 
               df_sample_data$days_wrt_delivery !=0))
length(which(df_sample_data$AP_or_PP == 1 & df_sample_data$data_available ==1 & 
               df_sample_data$days_wrt_delivery !=0))

length(which(df_sample_data$AP_or_PP == 0 & df_sample_data$data_available ==1 & 
               df_sample_data$days_wrt_delivery !=0 & df_sample_data$type == 'case'))
length(which(df_sample_data$AP_or_PP == 0 & df_sample_data$data_available ==1 & 
               df_sample_data$days_wrt_delivery !=0 & df_sample_data$type == 'control'))

length(which(df_sample_data$AP_or_PP == 1 & df_sample_data$data_available ==1 & 
               df_sample_data$days_wrt_delivery !=0 & df_sample_data$type == 'case'))
length(which(df_sample_data$AP_or_PP == 1 & df_sample_data$data_available ==1 & 
               df_sample_data$days_wrt_delivery !=0 & df_sample_data$type == 'control'))


