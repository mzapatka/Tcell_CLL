#with cut-off version

scDataMerged_harmony = NormalizeData(scDataMerged_harmony)
scDataMerged_harmony = ScaleData(scDataMerged_harmony)
#multimarkers
#exhaustion markers
#TOX, TIGIT, CTLA4, TNFRSF9 
gene.set <- c("TOX","PDCD1","CXCL13","LAG3","TIGIT","CTLA4","TNFRSF9", "EOMES", "ENTPD1", "HAVCR2", "CD38", "CD244A", "CD101")
scDataMerged_harmony <- AddModuleScore(
object = scDataMerged_harmony,
features = list(gene.set),
name = 'Exhaustion'
)
plot <- FeaturePlot(scDataMerged_harmony, 
					cols = c("gray", "red"),
          reduction = "umap", 
          features = c("Exhaustion1"),pt.size = 1.5, order=T,slot="scale.data",min.cutoff=0.3) + theme_minimal()
plot
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/TPEX/umap_ex_sig_",paste(gene.set,collapse="_"),"_V2.pdf"),width=13,height=8, plot=plot)



gene.set <- c("LAG3","XCL1","CRTAM","IFNG","CCL4","PDCD1","DUSP4","CD8A","ZEB2","NR4A2","SLA","NKG7","TIGIT","CTSW","TNFRSF9","TOX","LYST","TNFSF4","CCL3","GZMB","RAB27A","PRF1","CD70","PLSCR1","CXCL13")
scDataMerged_harmony <- AddModuleScore(
object = scDataMerged_harmony,
features = list(gene.set),
name = 'TPEX1'
)
plot <- FeaturePlot(scDataMerged_harmony, 
					cols = c("gray", "red"),
          reduction = "umap", 
          features = c("TPEX11"),pt.size = 1.5, order=T,slot="scale.data",min.cutoff=0.3)  + theme_minimal() 
plot
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/umap_TPEX_sig_1_V2.pdf"),width=13,height=8, plot=plot)

gene.set <- c("GZMK","CCL4L1","ITM2C","CD74","CCL4","AOAH","CXCR4","DTHD1","CCL3L3","CCL3L1","CLDND1","CD44","SH2D1A","TRAT1","EOMES","CCL5","F2R","TC2N","FAM102A","PVRIG","PDE4DIP","CMC1","CRTAM","NEK7","EPHA1","WIPF1","MS4A1","CD84")
scDataMerged_harmony <- AddModuleScore(
object = scDataMerged_harmony,
features = list(gene.set),
name = 'TPEX2'
)
plot <- FeaturePlot(scDataMerged_harmony, 
					cols = c("gray", "red"),
          reduction = "umap", 
          features = c("TPEX21"),pt.size = 1.5, order=T,slot="scale.data",min.cutoff=0.5)  + theme_minimal()
plot
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/umap_TPEX_sig_2_V2.pdf"),width=13,height=8, plot=plot)


#scDataMerged_harmony = SCTransform(scDataMerged_harmony, verbose = TRUE) 
dp = DotPlot(object = scDataMerged_harmony, features =c("TPEX21")) + RotatedAxis()
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/dot_TPEX_sig_2_V2.pdf"),width=13,height=8, plot=dp)


vp = VlnPlot(object = scDataMerged_harmony, features =c("Exhaustion1"), group.by="anno_laura_2022") + coord_flip() + NoLegend()
vp
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/vln_EX_sig_1_V3.pdf"),width=6,height=8, plot=vp)


vp = VlnPlot(object = scDataMerged_harmony, features =c("TPEX11"), group.by="anno_laura_2022" )  + coord_flip() + NoLegend()
vp
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/vln_TPEX_sig_1_V3.pdf"),width=6,height=8, plot=vp)


vp = VlnPlot(object = scDataMerged_harmony, features =c("TPEX21"), group.by="anno_laura_2022" )  + coord_flip() + NoLegend()
vp 
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/vln_TPEX_sig_2_V3.pdf"),width=6,height=8, plot=vp)



################################################################################
########## 			VERSION FOR EXHAUSTION MARKER PLOTS
#multimarkers
#exhaustion markers
#TOX, TIGIT, CTLA4, TNFRSF9 
gene.set <- c("TOX","PDCD1","CXCL13","LAG3","TIGIT","CTLA4","TNFRSF9", "EOMES", "ENTPD1", "HAVCR2", "CD38", "CD244A", "CD101")
scDataMerged_harmony <- AddModuleScore(
object = scDataMerged_harmony,
features = list(gene.set),
name = 'Exhaustion'
)
plot <- FeaturePlot(scDataMerged_harmony, 
					cols = c("gray", "red"),
          reduction = "umap", 
          features = c("Exhaustion1"),pt.size = 1.5, order=T,slot="scale.data") + theme_minimal()
plot
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/TPEX/umap_ex_sig_",paste(gene.set,collapse="_"),".pdf"),width=13,height=8, plot=plot)



gene.set <- c("LAG3","XCL1","CRTAM","IFNG","CCL4","PDCD1","DUSP4","CD8A","ZEB2","NR4A2","SLA","NKG7","TIGIT","CTSW","TNFRSF9","TOX","LYST","TNFSF4","CCL3","GZMB","RAB27A","PRF1","CD70","PLSCR1","CXCL13")
scDataMerged_harmony <- AddModuleScore(
object = scDataMerged_harmony,
features = list(gene.set),
name = 'TPEX1'
)
plot <- FeaturePlot(scDataMerged_harmony, 
					cols = c("gray", "red"),
          reduction = "umap", 
          features = c("TPEX11"),pt.size = 1.5, order=T,slot="scale.data")  + theme_minimal() 
plot
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/umap_TPEX_sig_1.pdf"),width=13,height=8, plot=plot)

gene.set <- c("GZMK","CCL4L1","ITM2C","CD74","CCL4","AOAH","CXCR4","DTHD1","CCL3L3","CCL3L1","CLDND1","CD44","SH2D1A","TRAT1","EOMES","CCL5","F2R","TC2N","FAM102A","PVRIG","PDE4DIP","CMC1","CRTAM","NEK7","EPHA1","WIPF1","MS4A1","CD84")
scDataMerged_harmony <- AddModuleScore(
object = scDataMerged_harmony,
features = list(gene.set),
name = 'TPEX2'
)
plot <- FeaturePlot(scDataMerged_harmony, 
					cols = c("gray", "red"),
          reduction = "umap", 
          features = c("TPEX21"),pt.size = 1.5, order=T,slot="scale.data",min.cutoff=0.5)  + theme_minimal()
plot
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/merged_CLL/TPEX/umap_TPEX_sig_2.pdf"),width=13,height=8, plot=plot)

