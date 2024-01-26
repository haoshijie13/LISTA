dir.create("pathway_phyper_2")
library(Seurat)
library(ggplot2)
SeuObj<-readRDS("D0.rds")
DY1=subset(SeuObj,subset=sample=="DY1_D0")
#DY1=subset(DY1,subset=coor_x>240 & coor_x<260 & coor_y>90 & coor_y<110)
SeuObj<-DY1
pathway=read.table("pathway_phyper_1.list",header=T,sep="\t")
pathway=pathway[pathway$Pvalue<0.01,]
res=apply(as.matrix(pathway),1,function(x){strsplit(x[2],",")[[1]]})
lens=apply(as.matrix(pathway),1,function(x){length(strsplit(x[2],",")[[1]])})
res=res[lens>2]
SeuObj<-AddModuleScore(SeuObj,features=res)
colnames(SeuObj@meta.data)=c(colnames(SeuObj@meta.data)[1:31],gsub(pattern="[ ,\\/-]+",replacement="_",sub(" \\[.*","",names(res))))
apply(as.matrix(colnames(SeuObj@meta.data)[32:ncol(SeuObj@meta.data)]),1,function(x){SeuObj@meta.data[,x][SeuObj@meta.data[,x]<quantile(SeuObj@meta.data[,x],0.005)[[1]]]=quantile(SeuObj@meta.data[,x],0.005)[[1]];SeuObj@meta.data[,x][SeuObj@meta.data[,x]>quantile(SeuObj@meta.data[,x],0.999)[[1]]]=quantile(SeuObj@meta.data[,x],0.999)[[1]];p<-ggplot(SeuObj@meta.data,aes(x=coor_x,y=coor_y,z=distance,fill=SeuObj@meta.data[,x]))+geom_tile()+geom_contour(bins=2,color="white")+coord_fixed()+theme_void()+scale_fill_viridis_c(option="B",direction=-1);ggsave(filename=paste0("pathway_phyper_2/",x,".pdf"),plot=p)})
#res1=res
#names(res1)=gsub(pattern="[ ,\\/-]+",replacement="_",sub(" \\[.*","",names(res)))
#apply(as.matrix(names(res1)),1,function(x){dir.create(paste0("pathway_phyper/",x));apply(as.matrix(res1[[x]]),1,function(y){df=SeuObj@meta.data;df$gene=SeuObj@assays$RNA@data[y,];df$gene[df$gene > quantile(df$gene,0.995)[[1]]]=quantile(df$gene,0.995)[[1]];p<-ggplot(df,aes(x=coor_x,y=coor_y,fill=gene))+geom_tile()+coord_fixed()+theme_void()+scale_fill_viridis_c(option="B",direction=-1);ggsave(filename=paste0("pathway_phyper/",x,"/",y,".pdf"),plot=p);return("none")});return("none")})
