args<-commandArgs(T)
library(data.table)
load(args[1])
sc=readRDS("/hwfssz1/ST_SUPERCELLS/P18Z10200N0350/99.lyw/liver/08.download/Cell_data/GSE192742.CD45.rds")
#sc=subset(sc,subset=annotation_lyw!="Erythrocyte")
library(Seurat)
#exp_spatial=FetchData(SeuObj,vars=rownames(SeuObj@assays$RNA@counts),slot="counts")
exp_spatial=as.matrix(SeuObj@assays$RNA@counts)
exp_spatial=as.data.frame(exp_spatial)
coord_spatial=SeuObj@meta.data[,c("coor_x","coor_y")]
nUMI_spatial=SeuObj@meta.data[,"nCount_RNA"]
names(nUMI_spatial)=rownames(SeuObj@meta.data)

#exp_sc=FetchData(sc,vars=rownames(sc@assays$RNA@counts),slot="counts")
exp_sc=as.matrix(sc@assays$RNA@counts)
exp_sc=as.data.frame(exp_sc)
celltype_sc=sc$annotation
celltype_sc=as.factor(celltype_sc)
nUMI_sc=sc@meta.data[,"nCount_RNA"]
names(nUMI_sc)=rownames(sc@meta.data)

library(RCTD)
reference <- Reference(exp_sc, celltype_sc, nUMI_sc)
puck <- SpatialRNA(coord_spatial, exp_spatial, nUMI_spatial)
rm(sc)
rm(SeuObj)
gc()
myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')

save(myRCTD,file="myRCTD_20220216.Rdata")

#results <- myRCTD@results
#norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/')
#cell_type_names <- myRCTD@cell_type_info$info[[2]]
#spatialRNA <- myRCTD@spatialRNA
#resultsdir <- 'RCTD_LVEC'
#dir.create(resultsdir)
#plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
#plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
#plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
#plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
#                     results$results_df)
#plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
#
#
#save(myRCTD,file="myRCTD_LVEC1.Rdata")
