args<-commandArgs(T)
library(Seurat)
SeuObj<-readRDS("add_Time_sample_change_marker_rank_factor.rds")  # Load the seurat object you want to deal with.
sub=subset(SeuObj,subset=Time==args[1])  # Select the dataset part for further analysis. OPTIONAL.
sub$group=as.integer(sample(1:nrow(sub@meta.data)/9/50,nrow(sub@meta.data),replace=T))  # Generate a random group label for each bin spot.
#lr<-read.table("../28.ligand/mouse_lr_pair.txt",header=T)
lr<-read.table("manual.lr.list",header=T)  # Read the ligand receptor list in this reportoiry or input a subset you want to deal with. 
genes=unique(c(lr[,2],lr[,3]))
genes=genes[genes %in% rownames(sub@assays$RNA@counts)]

#grp=cut(1:nrow(sub@assays$RNA@data),breaks = 100,labels = 1:100)
#x=split(rownames(sub@assays$RNA@data),f = grp)

exp_mat=aggregate(as.data.frame(t(as.matrix(sub@assays$RNA@counts[genes,]))),by=list("rank"=sub$rank,"group"=sub$group),mean)  # Aggregate the gene expression based on layer and group layer. Each layer will be split into random 50 groups.

res=apply(as.matrix(lr[ lr$ligand_gene_symbol %in% colnames(exp_mat) & lr$receptor_gene_symbol %in% colnames(exp_mat),]),1,function(x){x=as.vector(x); exp_mat[,x[2]] * exp_mat[,x[3]] })  # Multiply the ligand gene expression by the corresponded receptor gene expression of each group as the interaction strength score.
colnames(res)=lr[ lr$ligand_gene_symbol %in% colnames(exp_mat) & lr$receptor_gene_symbol %in% colnames(exp_mat),"lr_pair"]
res=as.data.frame(res)
res$rank=exp_mat$rank
pvalues=apply(as.matrix(colnames(res)[1:(ncol(res)-1)]),1,function(x){test=kruskal.test( get(x) ~ rank, data = res);test$p.value})    # KW test
names(pvalues)=colnames(res)[1:(ncol(res)-1)]
test=aggregate(res,by=list("rank"=res$rank),mean)
test=as.data.frame(test)

pdf(paste0("manual_lxr.",args[1],".pdf"),height=length(na.omit(pvalues[pvalues<0.05]))/5.5,width=5)
pheatmap::pheatmap(t(test[,names(na.omit(pvalues[pvalues<0.05]))]),scale="row",cluster_cols=F,cluster_rows=T)
dev.off()
#df=as.data.frame(cbind(t(test[,names(na.omit(pvalues[pvalues<0.05]))]),na.omit(pvalues[pvalues<0.05])))
df=as.data.frame(cbind(t(test[,names(pvalues)]),pvalues))
colnames(df)=c(paste0("layer_",1:9),"Pvalue")
write.table(df,file=paste0("manual_lxr.",args[1],".table.xls"),sep="\t",quote=F,row.names=T,col.names=T)
