
# find which genes in NT are different than all other samples (pairwise comparisons)

res_NT_BAT<-subset(results(dds, contrast=c("cell", "NT", "BAT"), alpha=0.01, lfcThreshold=1), padj < 0.01)
res_NT_BMDM<-subset(results(dds, contrast=c("cell", "NT", "BMDM"), alpha=0.01, lfcThreshold=1), padj < 0.01)
res_NT_Cortex<-subset(results(dds, contrast=c("cell", "NT", "Cortex"), alpha=0.01, lfcThreshold=1), padj < 0.01)
res_NT_MEF<-subset(results(dds, contrast=c("cell", "NT", "MEF"), alpha=0.01, lfcThreshold=1), padj < 0.01)

NT_BAT_down<-rownames(head(res_NT_BAT[ order(res_NT_BAT$log2FoldChange), ], floor(length(rownames(res_NT_BAT))/ 50)))
NT_BMDM_down<-rownames(head(res_NT_BMDM[ order(res_NT_BMDM$log2FoldChange), ], floor(length(rownames(res_NT_BMDM))/ 50)))
NT_Cortex_down<-rownames(head(res_NT_Cortex[ order(res_NT_Cortex$log2FoldChange), ], floor(length(rownames(res_NT_Cortex))/50)))
NT_MEF_down<-rownames(head(res_NT_MEF[ order(res_NT_MEF$log2FoldChange), ], floor(length(rownames(res_NT_MEF))/50)))

# find common downregulated genes

down<-Reduce(intersect, list(NT_BAT_down, NT_BMDM_down, NT_Cortex_down, NT_MEF_down))

NT_BAT_up<-rownames(head(res_NT_BAT[ order(res_NT_BAT$log2FoldChange, decreasing=TRUE), ], floor(length(rownames(res_NT_BAT))/ 50)))
NT_BMDM_up<-rownames(head(res_NT_BMDM[ order(res_NT_BMDM$log2FoldChange, decreasing=TRUE), ], floor(length(rownames(res_NT_BMDM))/ 50)))
NT_Cortex_up<-rownames(head(res_NT_Cortex[ order(res_NT_Cortex$log2FoldChange, decreasing=TRUE), ], floor(length(rownames(res_NT_Cortex))/50)))
NT_MEF_up<-rownames(head(res_NT_MEF[ order(res_NT_MEF$log2FoldChange, decreasing=TRUE), ], floor(length(rownames(res_NT_MEF))/50)))



# find common upregulated genes

up<-Reduce(intersect, list(NT_BAT_up, NT_BMDM_up, NT_Cortex_up, NT_MEF_up))

length(down)
length(up)

# make heatmap

pdf("gene_heatmap_DOWN.pdf", height=20)
matDown<-assay(rld)[down,]
matDown<-matDown - rowMeans(matDown)
dfDown<-as.data.frame(colData(rld)[ "cell" ])
pheatmap(matDown, annotation_col=dfDown)

dev.off()


pdf("gene_heatmap_UP.pdf", height=20)
matUp<-assay(rld)[up,]
matUp<-matUp - rowMeans(matUp)
dfUp<-as.data.frame(colData(rld)[ "cell" ])
pheatmap(matUp, annotation_col=dfUp)

dev.off()

pdf("gene_heatmap_ALL.pdf", height=50)
matAll<-assay(rld)
matAll<-matAll - rowMeans(matAll)
dfAll<-as.data.frame(colData(rld)[ "cell" ])
pheatmap(matAll, annotation_col=dfAll)

dev.off()

