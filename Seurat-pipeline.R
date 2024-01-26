library("Seurat")
library("ggplot2")
library("reshape2")
library("tidyr")
library(data.table)
library(Matrix)
# library("future")
# plan("multiprocess", workers = 30)
options(future.globals.maxSize=100000000000)
args=commandArgs(T)
base=basename(args[1])
dir=dirname(args[1])
prefix=sub(".txt","",base)

read<-function(mat,bin){
        data<-fread(mat,header=T)
        data$cellID<-paste(as.character(round(data$x/bin, digits = 0)),"_",as.character(round(data$y/bin, digits = 0)),sep="")
        gene=unique(data$geneID)
        cell=unique(data$cellID)
        gene_idx=c(1:length(gene))
        cell_idx=c(1:length(cell))
        names(gene_idx)=gene
        names(cell_idx)=cell
        print(head(gene_idx[data$geneID]))
        data=as.data.frame(data)
        mat=sparseMatrix(i=gene_idx[data$geneID],j=cell_idx[data$cellID],x=data[,4])
        rownames(mat)=gene
        colnames(mat)=cell
        return(mat)
}

analyze <- function(mat,bin){
        data=mat
        data[is.na(data)]=0
        print("ok")
        SeuObj<-CreateSeuratObject(counts = data, names.delim = "-", project = "SeuObj")
        SeuObj[["percent.mt"]] <- PercentageFeatureSet(SeuObj, pattern = "^mt-|^MT-|^Mt-")
        #save(SeuObj,file=paste(dir,"/bin",bin,".Rdata",sep=""))
        #q()
        SeuObj <- NormalizeData(SeuObj, normalization.method = "LogNormalize", scale.factor = 10000)
        SeuObj <- FindVariableFeatures(SeuObj, selection.method = "vst", nfeatures = 2000)
        all.genes <- rownames(SeuObj)
        SeuObj <- ScaleData(SeuObj, features = all.genes, vars.to.regress = "nCount_RNA")
#        SeuObj <- SCTransform(SeuObj, vars.to.regress = "percent.mt", verbose = FALSE)
        SeuObj <- RunPCA(SeuObj, verbose = FALSE)
        SeuObj <- RunUMAP(SeuObj, dims = 1:15, verbose = FALSE)
        SeuObj <- FindNeighbors(SeuObj, dims = 1:15, verbose = FALSE)
        SeuObj <- FindClusters(SeuObj, verbose = FALSE)
        SeuObj@meta.data$coor_x=sub(rownames(SeuObj@meta.data),pattern = "_.*",replacement = "")
        SeuObj@meta.data$coor_y=sub(rownames(SeuObj@meta.data),pattern = ".*_",replacement = "")
        SeuObj@meta.data$coor_x=sub(SeuObj@meta.data$coor_x,pattern = "X",replacement = "")
        SeuObj@meta.data$coor_x=as.integer(SeuObj@meta.data$coor_x)
        SeuObj@meta.data$coor_y=as.integer(SeuObj@meta.data$coor_y)
        AllMG <<- FindAllMarkers(SeuObj)
        SeuObj <<- SeuObj
#       if(bin==100){
#               print(mean(SeuObj@meta.data$nFeature_RNA))
#               print(mean(SeuObj@meta.data$nCount_RNA))
#               count = as.data.frame(SeuObj@assays$RNA@counts)
#               count_1=count>0
#               count_1=as.data.frame(count_1)
#               print(mean(rowSums(as.data.frame(count_1)))/ncol(count_1))
#       }
    filename=paste(dir,"/bin",bin,".MG",sep="")
    write.table(AllMG,file=filename)
    save(SeuObj,file=paste(dir,"/bin",bin,".Rdata",sep=""))
}

bins=c(100,50)
for(bin in bins){
        mat<-read(args[1],bin)
        analyze(mat,bin)
}
