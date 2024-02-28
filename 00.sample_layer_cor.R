library(Seurat)
library(dplyr)
SeuObj<-readRDS("add_Time_sample_change_marker_rank_factor.rds") # load rds file including all sections 
grp=cut(1:nrow(SeuObj@assays$RNA@data),breaks = 100,labels = 1:100)
x=split(rownames(SeuObj@assays$RNA@data),f = grp)
#result=lapply(x,function(x){exp=FetchData(SeuObj,vars = x);aggregate(exp,by=list("layer"=paste0(SeuObj$Time,"_",SeuObj$rank)),mean)})
#test<-bind_cols(result, .name_repair="unique")
#saveRDS(test,file="layer_exp.rds")

result=lapply(x,function(x){exp=FetchData(SeuObj,vars = x);aggregate(exp,by=list("layer"=paste0(SeuObj$sample,"_",SeuObj$rank)),mean)})
test<-bind_cols(result, .name_repair="unique")
saveRDS(test,file="sample_layer_exp.rds")

#result=lapply(x,function(x){exp=FetchData(SeuObj,vars = x);aggregate(exp,by=list("layer"=SeuObj$sample),mean)})
#test<-bind_cols(result, .name_repair="unique")
#saveRDS(test,file="sample_exp.rds")

#df<-readRDS("sample_layer_exp.rds")
rownames(df)=df[,1]
df=df[,-1]
df=t(df)
res=apply(df,2,as.numeric)
res=as.data.frame(res)
rownames(res)=rownames(df)
res=res[complete.cases(res),]
cor_res=cor(res)
pdf("sample_layer_cor.pdf",width=30,height=30)
pheatmap::pheatmap(cor_res[order(sub("_.*_","_",sub("[^_]*_","",rownames(cor_res)))),][,order(sub("_.*_","_",sub("[^_]*_","",colnames(cor_res))))][c(91:135,37:90,136:306),c(91:135,37:90,136:306)],cluster_rows=F,cluster_cols=F)
dev.off()
write.table(cor_res[order(sub("_.*_","_",sub("[^_]*_","",rownames(cor_res)))),][,order(sub("_.*_","_",sub("[^_]*_","",colnames(cor_res))))][c(91:135,37:90,136:306),c(91:135,37:90,136:306)],file="sample_layer_cor.xls",row.names=T,col.names=T,sep="\t")
