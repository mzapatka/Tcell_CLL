
cellchat.NL <- readRDS("cellchat/CLL_cellchatliana_L2_CD4_8_A.rds")
cellchat.LS <- readRDS("cellchat/rLN_cellchatliana_L2_CD4_8_A.rds")
UpdateIdents=FALSE
LiftIdents=FALSE

if (LiftIdents){
	group.new = unique(c("B",as.character(cellchat.NL@idents),as.character(cellchat.LS@idents)))
	#group.new = levels(cellchat.NL@idents)
	cellchat.LS <- liftCellChat(cellchat.LS, group.new)
	cellchat.NL <- liftCellChat(cellchat.NL, group.new)
}else if (UpdateIdents) {
	cellchat.LS@idents = factor(as.character(cellchat.LS@idents), levels=levels(cellchat.NL@idents))
}


object.list <- list(CLL = cellchat.NL, rLN = cellchat.LS)
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))#, cell.prefix = TRUE)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

pdf(paste("circos_D_interactions.pdf",sep=""),width=8,height=8)
#png(cellchat/compare_D_interactions.png,width=1300,height=750,res=150)
par(mfrow = c(1,1))#, xpd=TRUE)
#cellchat@idents$joint = factor(cellchat@idents$joint,levels=c("CD4 TCM","Tregs","CD8 TEM","CD4 TN","CD8 TN","TFH","CD4 TEM","B"))
netVisual_diffInteraction(cellchat, weight.scale = T, color.edge = c( "#2166ac", "#b2182b")) 
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", color.edge = c( "#2166ac", "#b2182b"))
#netVisual_diffInteraction(cellchat, weight.scale = T,source.use=rev(seq(1,10,1)),targets.use	=rev(seq(1,10,1)))
dev.off()


#diff interactions
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf(paste("/cellchat/compare_n_interactions.pdf",sep=""),width=10,height=10)#,res=150)
gg1 + gg2
dev.off()

#plot N interactions
pdf(paste("cellchat/circos_n_interactions.pdf",sep=""),width=10,height=10)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

#IN/OUT interaction HM
gg1 <- netVisual_heatmap(cellchat, comparison=c(2,1))#), color.heatmap = c( "#b2182b", "#2166ac"))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, comparison=c(2,1), measure = "weight")# , color.heatmap = c( "#b2182b", "#2166ac"))
#> Do heatmap based on a merged object
pdf(paste("cellchat/n_interaction_HM.pdf",sep=""),width=9,height=5)
gg1 + gg2
dev.off()


#pairwise
#levelsuse = c(2:length(levels(cellchat@idents$joint)))
labelsuseContrast = c("CD4 TCM","CD4 TEM","CD8 TEM","CD4 TN","CD8 TN","TFH","Tregs")
levelsuseContrast = c(which(levels(cellchat@idents$joint) %in% labelsuseContrast))
levelsuseCase = c(which(levels(cellchat@idents$joint) %in% "B"))

pdf(paste("cellchat/pairwise_dot_diff_interactions.pdf",sep=""),width=10,height=10)#,res=150)
pairwise_dot_diff_interactions <- netVisual_bubble(cellchat, sources.use = levelsuseCase, targets.use = levelsuseContrast,  comparison = c(1, 2), angle.x = 45)
pairwise_dot_diff_interactions
dev.off()
write.table(pairwise_dot_diff_interactions$data,paste("cellchat/pairwise_dot_diff_interactions.tsv",sep=""),sep="\t")

#pairwise rev
pdf(paste("cellchat/pairwise_dot_diff_interactions_rev.pdf",sep=""),width=10,height=10)#,res=150)
pairwise_dot_diff_interactions_rev<-netVisual_bubble(cellchat, sources.use = levelsuseContrast, targets.use = levelsuseCase,  comparison = c(1, 2), angle.x = 45)
pairwise_dot_diff_interactions_rev
dev.off()
write.table(pairwise_dot_diff_interactions_rev$data,paste("cellchat/pairwise_dot_diff_interactions_rev.tsv",sep=""),sep="\t")



