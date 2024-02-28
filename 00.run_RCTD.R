args<-commandArgs(T)
library(data.table)
load(args[1])   # Load the seurat object of your sptaial data which named by SeuObj
sc=readRDS("GSE192742.CD45.rds") # Load the annotated single cell data 
#sc=subset(sc,subset=annotation_lyw!="Erythrocyte")
library(Seurat) 
#exp_spatial=FetchData(SeuObj,vars=rownames(SeuObj@assays$RNA@counts),slot="counts")
exp_spatial=as.matrix(SeuObj@assays$RNA@counts) 
exp_spatial=as.data.frame(exp_spatial) 
coord_spatial=SeuObj@meta.data[,c("coor_x","coor_y")] # the coordinations of each spot were stored in variables "coor_x" and "coor_y"
nUMI_spatial=SeuObj@meta.data[,"nCount_RNA"]
names(nUMI_spatial)=rownames(SeuObj@meta.data)

#exp_sc=FetchData(sc,vars=rownames(sc@assays$RNA@counts),slot="counts")
exp_sc=as.matrix(sc@assays$RNA@counts)
exp_sc=as.data.frame(exp_sc)
celltype_sc=sc$annotation  # Specify your single cell annotation
celltype_sc=as.factor(celltype_sc)
nUMI_sc=sc@meta.data[,"nCount_RNA"]
names(nUMI_sc)=rownames(sc@meta.data)

library(RCTD)
# Create RCTD object
reference <- Reference(exp_sc, celltype_sc, nUMI_sc)
puck <- SpatialRNA(coord_spatial, exp_spatial, nUMI_spatial)
# Clean the environment
rm(sc)
rm(SeuObj)
gc()
# Run RCTD analysis
myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')

save(myRCTD,file="myRCTD_20220216.Rdata")
