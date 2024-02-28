library(DESeq2)
library(ggplot2)
library(pheatmap)
library(data.table)


count_final <- read.delim("TBLXR1_count.txt")
head(count_final)

samplenames <- colnames(count_final)
group <- c("NC", "NC","NC","NC","KD","KD","KD","KD")
count_final <- as.matrix(count_final)
table.all <- data.frame(name = samplenames, 
                        condition=group)
dds.all <- DESeqDataSetFromMatrix(floor(count_final), colData=table.all, design= ~ condition)
dds.all <- dds.all[ rowSums(counts(dds.all)) > 1, ]
deg = results(dds.all, contrast = c("condition","KD","NC"))
deg = deg[deg$pvalue < 0.05,]
deg_up = deg[deg$log2FoldChange >  log(1.5, 2),]
deg_down = deg[deg$log2FoldChange <  -log(1.5, 2),]
write.csv(deg_up, 'KD_tbl1xr1.up_FC1.5.csv', row.names = F)
write.csv(deg_down, 'KD_tbl1xr1.down_FC1.5.csv', row.names = F)

deg$change = 'No'
deg$change[match(deg_up$SYMBOL, deg$SYMBOL)] = 'Up'
deg$change[match(deg_down$SYMBOL, deg$SYMBOL)] = 'Down'

deg_down_labeling = unique(deg_down$SYMBOL)
gene = c('Ccnd1','Tbl1xr1','Axin2','Acadm','Crot','Lgr5')
deg_down_labeling_sel = deg_down_labeling[match(gene, deg_down_labeling$SYMBOL),]

p <- ggplot(data=deg, aes(x=log2FoldChange, y=log10pvalue, col = change, label = labeling)) + geom_point()+ theme_minimal()+   
scale_color_manual(values=c("blue", "black", "red"))+ geom_text( color="black")

ggsave('volcano_plot.pdf',p)
