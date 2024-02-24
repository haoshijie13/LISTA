pathway=read.table("all.kegg.pathway.list",sep="\t")            # Input gene lists of pathways.
pathways=apply(as.matrix(pathway$V2),1,function(x){strsplit(x,split=",")[[1]]})
names(pathways)=pathway$V1

zonated_gene=read.table("zonated.gene_Halpern.list")            # Input zonation gene list
zonated_gene=zonated_gene$V1 
exp=readRDS("layer_exp_change.rds")                             # Input layer averaged gene expression
zonated_gene=zonated_gene[zonated_gene %in% rownames(exp)]
library(dplyr)
pvalues=lapply(pathways,function(x){1-phyper(length(x[x %in% zonated_gene]),length(zonated_gene),nrow(exp),length(x))})    # Hypergeometric test of genes of each pathway against zonation gene
df=t(as.data.frame(pvalues))
rownames(df)=names(pvalues)
#write.table(df[df[,1]<0.05,],sep="\t",quote=F,file="pathway_phyper.list")
pathw_gene=lapply(pathways,function(x){paste0(x[x %in% zonated_gene],collapse=",")})
df1=t(as.data.frame(pathw_gene))
rownames(df1)=names(pathw_gene)
df2=merge(df,df1,by=0)
rownames(df2)=df2$Row.names
df2=df2[,-1]
colnames(df2)=c("Pvalue","genes")
#write.table(df2[df2[,1]<0.05,],sep="\t",quote=F,file="pathway_phyper.list")
write.table(df2[df2[,1]<0.05,],sep="\t",quote=F,file="pathway_phyper_1.list")
res=apply(as.matrix(df2[df2$Pvalue<0.05,]),1,function(x){genes=strsplit(x[2],split=",")[[1]];if(length(genes)>1){colSums(exp[genes,19:27])}else{exp[genes,19:27]}})
res1=as.data.frame(t(res))

pathway_gene=lapply(pathways,function(x){paste0(x,collapse=",")})
tab=cbind(df2[rownames(res1),],res1,as.character(pathw_gene[rownames(res1)]),as.character(pathway_gene[rownames(res1)]))
colnames(tab)=c("Pvalue","zonated genes",colnames(tab)[3:11],"zonated genes","all genes")
tab=tab[,-2]
tab$number_zg=apply(as.matrix(tab$`zonated gene`),1,function(x){length(strsplit(x,",")[[1]])})
tab$number_gene=apply(as.matrix(tab$`all genes`),1,function(x){length(strsplit(x,",")[[1]])})
tab$qvalue=qvalue::qvalue(tab$Pvalue, lambda = seq(0, max(tab$Pvalue), 0.05))$qvalue
tab=tab[,c(1,15,13,14,2:10,11,12)]
write.table(tab,row.names=T,file="pathway_phyper_table_halpern.tsv",sep="\t",quote=F)
