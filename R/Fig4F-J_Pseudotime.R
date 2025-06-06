#General Steps:
#Setup and Data Loading: Loads necessary libraries and reads in a pre-processed SingleCellExperiment object (scDataMerged_harmony).
#Pseudotime Calculation: Uses destiny to calculate pseudotime (diffusion maps and DPT).
#Differential Expression Analysis: Uses tradeSeq to identify genes changing along the pseudotime trajectory.
#Visualization: Creates various plots, including diffusion maps, heatmaps, and pseudotime plots, to visualize the results.
#Clustering and Heatmap Refinement: Performs k-means clustering on the differentially expressed genes and refines the heatmaps based on the clustering results.


library(patchwork)
library(gg3D)
#library(plotly)
library(ggthemes )
require(destiny)
library(SingleCellExperiment)
library(gridExtra)
library(rgl)
library(RColorBrewer)
library(BiocParallel)
library(ComplexHeatmap)


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

projectdir="CLL_TCLmouse"
mergeID="merged_CLL"
scDataMerged_harmony=readRDS(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/seurat_harmony_umap_obj_formonocle3.rds"))


##################################################
#CD4
CD4clusts = c(0,1,4,5,6,10)
scDataMerged.CD4 = subset(scDataMerged_harmony,subset = seurat_clusters %in% CD4clusts)
#scDataMerged.CD4 <- SCTransform(scDataMerged.CD4, verbose = TRUE, variable.features.n = dim(scDataMerged.CD4@assays$RNA)[1])


expr = as.matrix(scDataMerged.CD4@assays$integrated@data) # the expression matrix
types = as.character(scDataMerged.CD4@meta.data$seurat_clusters) # a character list of the types. Samples from the same type should have the same name.
main_types = as.character(unique(types)) # a character list of the main types. 

#expr_SE <- SingleCellExperiment(list(counts=expr))
#expr_SE <- logcounts(expr_SE)
#expr_SE@assays@data$logcounts = expr_SE@assays@data$counts	#disable lognorm on Integrated data

#expr_SE = as(expr_SE, "SingleCellExperiment")

colnames(expr) <- types

sce<-SingleCellExperiment(list(logcounts=expr),metadata=list(cluster=types))
dm <- DiffusionMap(sce,n_pcs = 50,verbose=TRUE)
dpt <- DPT(dm, tips = c(3986,2286,4993))

#plot 3D
color3d=brewer.pal(8, "Dark2")[1:6]
color3dlevels=data.frame(row.names=levels(factor(types)),color3d)
##preserve the umap color
##color3dlevelsmapped = color3dlevels[factor(types),]
#brand new coloring
colorhighcontrast = rev(brewer.pal(6, "Dark2"))


names(color3dlevelsmapped) = factor(types)
plot(dm, col=color3dlevelsmapped,pch=20,draw_legend=TRUE,interactive=TRUE)
legend3d("topright", legend = paste('Cluster', rownames(color3dlevels)), pch = 16, col = color3dlevels$color3d, cex=1, inset=c(0.02))

plot.DiffusionMap(dm, dims=2:4,col=color3dlevelsmapped,pch=20,draw_legend=TRUE,interactive=TRUE)



#plot 2D
#regular ggplot
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  Timepoint = types)
tmp$size <- 0.5
ggcolors = gg_color_hue(13)

#tmp$Timepoint = factor(tmp$Timepoint,levels=c("5","0","6","1","4","10"))
tmp$rand = runif(dim(tmp)[1], 1, dim(tmp)[1])

diffplot1 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint,size=size)) +
    geom_point(size=1) + scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic()+
    scale_color_manual(values=colorhighcontrast[c(1,2,3,6,4,5)]) 
    #scale_color_manual(values=ggcolors[c(1,2,11,5,6,7)]) +
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/destiny_CD4_diff_map12_C2.pdf"),width=13,height=8, plot=diffplot1)


diffplot1 <- ggplot(tmp %>% arrange(desc(rand)), aes(x = DC1, y = DC2, colour = Timepoint,size=size)) +
    geom_point(size=1) + scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic() +
    xlim(-0.08,0.02)+
    ylim(-0.07,0.07)+
    scale_color_manual(values=colorhighcontrast[c(1,2,3,6,4,5)]) 
    #scale_color_manual(values=ggcolors[c(1,2,11,5,6,7)]) +
diffplot1    
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/destiny_CD4_diff_map12zoom_C2.pdf"),width=13,height=8, plot=diffplot1)

#super zoom
diffplot1z <- ggplot(tmp %>% arrange(desc(rand)), aes(x = DC1, y = DC2, colour = Timepoint,size=size)) +
    geom_point(size=2) + scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic() +
    xlim(-0.015,0.015)+
    ylim(-0.015,0.015) +
    scale_color_manual(values=colorhighcontrast[c(1,2,3,6,4,5)]) 
    #scale_color_manual(values=ggcolors[c(1,2,11,5,6,7)]) 
diffplot1z   
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/destiny_CD4_diff_map12zoomB_C2.pdf"),width=13,height=8, plot=diffplot1z)


diffplot2 <- ggplot(tmp, aes(x = DC2, y = DC3, colour = Timepoint)) +
    geom_point() + scale_color_tableau() + 
    xlab("Diffusion component 2") + 
    ylab("Diffusion component 3") +
    theme_classic()


##################################################
#CD8

CD8clusts = c(2,3,7,8)
#CD8clusts = c(2,7,8)
scDataMerged.CD8 = subset(scDataMerged_harmony,subset = seurat_clusters %in% CD8clusts)
#scDataMerged.CD8 <- SCTransform(scDataMerged.CD8, verbose = TRUE, variable.features.n = dim(scDataMerged.CD8@assays$RNA)[1])


expr = as.matrix(scDataMerged.CD8@assays$integrated@data) # the expression matrix
types = as.character(scDataMerged.CD8@meta.data$seurat_clusters) # a character list of the types. Samples from the same type should have the same name.
main_types = as.character(unique(types)) # a character list of the main types. 

#expr_SE <- SingleCellExperiment(list(counts=expr))
#expr_SE <- logcounts(expr_SE)
#expr_SE@assays@data$logcounts = expr_SE@assays@data$counts	#disable lognorm on Integrated data

#expr_SE = as(expr_SE, "SingleCellExperiment")
phenoDF = AnnotatedDataFrame(data=scDataMerged.CD8@meta.data)
expr_ES =ExpressionSet(assayData=expr,phenoData=phenoDF)

colnames(expr) <- types

#dm <- DiffusionMap(t(expr),n_pcs = 50,verbose=TRUE)
dm <- DiffusionMap(expr_ES,n_pcs = 18,verbose=TRUE)
#dptCD8 <- DPT(dm)
#dptCD8 <- DPT(dm, tips = c(1911,2536,1431))	#4 clust, TEM1
#dptCD8 <- DPT(dm, tips = c(1911,2267,1431))	#4 clust, TEM2
dptCD8 <- DPT(dm, tips = c(1259,1672,245))	#3 clust, TEM1
plot(dptCD8,   col_by = 'Branch')
plot(dptCD8,   col_by = 'dpt')

#plot(dptCD8,   col_by = 'seurat_clusters')

#exhaustion
scDataMerged.CD8$DC1 = dm$DC1
scDataMerged.CD8$DC2 = dm$DC2
scDataMerged.CD8$DC3 = dm$DC3
plot(scDataMerged.CD8$DC1,scDataMerged.CD8$DC2,pch=20,col=scDataMerged.CD8$Exhaustion1)
require(ggplot2)
qplot(scDataMerged.CD8$DC1,scDataMerged.CD8$DC2,  colour = scDataMerged.CD8$Exhaustion1)

#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
library( colorRamps )
myPalette <- green2red
green2red(n)
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(scDataMerged.CD8$Exhaustion1), max(scDataMerged.CD8$Exhaustion1)))
metadataCD8=scDataMerged.CD8@meta.data
ggplot(metadataCD8, aes(x = DC1, y = DC2, colour = Exhaustion1)) +
  geom_point() +
  scale_size_area() + sc + theme_bw()

theta=110
phi=260
dplot3d <- ggplot(metadataCD8, aes(x = DC1, y = DC2, z = DC3, colour = Exhaustion1)) +
     scale_color_tableau() + 
    xlab("ignore axis 1") + 
    ylab("ignore axis 2")  +
  axes_3D(theta=theta,phi=phi) +
  labs_3D(labs=c("DC1", "DC2", "DC3"),theta=theta,phi=phi)+#vjust=c(-15,0,55),hjust=c(0,0,5) ) +
  stat_3D(theta=theta, phi=phi,size=0.75) + sc +
  theme_bw()
print(dplot3d)
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/pseudotime_destiny/destiny_CD8_diff_map3D_ex_score.pdf"),width=13,height=8, plot=dplot3d)


dm$ex_score = 1
plot(dm,col_by = 'dpr')
plot(dpt, root=3, 2:3, col_by = 'dpr')

#plot to pdf
ggcolors = gg_color_hue(13)
CD8ClusterIdx = c(3,4,8,9)
ggcolorCD8 = data.frame(row.names=c("2","3","7","8"),col=ggcolors[CD8ClusterIdx])


plot.DiffusionMap(dm, dims=1:3,pch=20,draw_legend=TRUE,interactive=FALSE,angle = 45)
#  You have 15189 genes. Consider passing e.g. n_pcs = 50 to speed up computation.




tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  Timepoint =types)
dplot <- ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
    geom_point() + scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2") +
    theme_classic()
dplot
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/pseudotime_destiny/destiny_CD8_diff_map.pdf"),width=13,height=8, plot=dplot)


#############
##CD8 3D plot
ggcolors = gg_color_hue(13)
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
                  Timepoint = types)
for(i in c(1:15)){
	#theta=310 #working but flipped
	theta=110
	#theta=-i * 10
	print(theta)
	#theta=i * 10
	phi=260
	#phi=i * 30 + 100
	print(phi)
	#https://karobben.github.io/2020/06/19/R/GG3D/#1-%E5%AE%89%E8%A3%85

	#theta=-80
	#phi=50
	dplot3d <- ggplot(tmp, aes(x = DC1, y = DC2, z = DC3, colour = Timepoint)) +
	     scale_color_tableau() + 
	    xlab("ignore axis 1") + 
	    ylab("ignore axis 2")  +
	  axes_3D(theta=theta,phi=phi) +
	  labs_3D(labs=c("DC1", "DC2", "DC3"),theta=theta,phi=phi)+#vjust=c(-15,0,55),hjust=c(0,0,5) ) +
	  stat_3D(theta=theta, phi=phi,size=0.75) +
    scale_color_manual(values=ggcolors[CD8ClusterIdx])  +
	  theme_bw()
	print(dplot3d)
}
ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/pseudotime_destiny/destiny_CD8_diff_map3D.pdf"),width=13,height=8, plot=dplot3d)


tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  Timepoint = types)
fig <- plot_ly(tmp, x = ~DC1, y = ~DC2, z = ~DC3, color = ~Timepoint)

fig


#plot 3d 2-4
plot.DiffusionMap(dm, dims=2:4,col=color3dlevelsmapped,pch=20,draw_legend=TRUE,interactive=TRUE)



######################################
#optimization

#  You have 15189 genes. Consider passing e.g. n_pcs = 50 to speed up computation.
plot(dm)

#look at other sigma values
sigmas <- find_sigmas(sce, verbose = FALSE)
optimal_sigma(sigmas)
for (sigma in c(10, 25, optimal_sigma(sigmas), 100))
	plot(DiffusionMap(sce, sigma,n_pcs = 50), 1:2,
	main = substitute(sigma == s, list(s = round(sigma,2))),
	col=types, draw.legend = FALSE)
	
for (sigma in c(10, 25, optimal_sigma(sigmas), 100))
	plot(DiffusionMap(sce, sigma,n_pcs = 50), 1:2,
	main = substitute(sigma == s, list(s = round(sigma,2))),
	col.by = 'num.cells', draw.legend = FALSE)


color3d=brewer.pal(8, "Dark2")[1:4]
color3dlevels=data.frame(row.names=levels(factor(types)),color3d)
color3dlevelsmapped = color3dlevels[factor(types),]
names(color3dlevelsmapped) = factor(types)
plot(dm, col=color3dlevelsmapped,pch=20,draw_legend=TRUE,interactive=TRUE)
legend3d("topright", legend = paste('Cluster', rownames(color3dlevels)), pch = 16, col = color3dlevels$color3d, cex=1, inset=c(0.02))

###############################

theta=1
#theta=i * 10
phi=1
#phi=i * 10
#https://karobben.github.io/2020/06/19/R/GG3D/#1-%E5%AE%89%E8%A3%85

dplot3d <- ggplot(tmp, aes(x = DC1, y = DC2, z = DC3, colour = Timepoint)) +
     scale_color_tableau() + 
    xlab("Diffusion component 1") + 
    ylab("Diffusion component 2")  +
  axes_3D(theta=theta,phi=phi) +
  stat_3D(theta=theta, phi=phi) +
  axis_labs_3D(theta=theta, phi=phi, size=1)
print(dplot3d)

#################
# 	dpt analysis
#################
#CD4
library(ggbeeswarm)
scDataMerged.CD4$pseudotime_dpt <- rank(dpt$dpt) 
scDataMerged.CD4$pseudotime_dpt = max(scDataMerged.CD4$pseudotime_dpt) - scDataMerged.CD4$pseudotime_dpt
DPT_by_celltype <- ggplot(scDataMerged.CD4@meta.data, 
       aes(x = pseudotime_dpt, 
           y = seurat_clusters, colour = seurat_clusters)) +
    scale_color_tableau() + theme_classic() +
    xlab("Diffusion map pseudotime (dpt)") +
    geom_quasirandom(groupOnX = FALSE) +
    ylab("Timepoint") +
    ggtitle("Cells ordered by diffusion map pseudotime")
DPT_by_celltype

ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/pseudotime_destiny/destiny_CD4_diff_steps.pdf"),width=13,height=8, plot=DPT_by_celltype)


#CD8
library(ggbeeswarm)
scDataMerged.CD8$pseudotime_dpt <- rank(dptCD8$dpt) 
scDataMerged.CD8$pseudotime_dpt = max(scDataMerged.CD8$pseudotime_dpt) - scDataMerged.CD8$pseudotime_dpt
DPT_by_celltype <- ggplot(scDataMerged.CD8@meta.data, 
       aes(x = pseudotime_dpt, 
           y = seurat_clusters, colour = seurat_clusters)) +
    scale_color_tableau() + theme_classic() +
    xlab("Diffusion map pseudotime (dpt)") +
    geom_quasirandom(groupOnX = FALSE) +
    ylab("Timepoint") +
    ggtitle("Cells ordered by diffusion map pseudotime")
DPT_by_celltype

ggplot2::ggsave(filename = paste(sep="",projectdir,"/analysis/markers/",mergeID,"/pseudotime_destiny/destiny_CD8_diff_steps.pdf"),width=13,height=8, plot=DPT_by_celltype)

#################
# 	dpt DEG analysis : simple
#################
library(tidyverse)

order_tip_1 <- order(dpt[tips(dpt)[[1]], ])

heatmaps <-
    na.exclude(unique(dpt$Branch)) %>%
    set_names() %>%
    map(function(branch) {
        branch_subset <- dpt$Branch == branch & !is.na(dpt$Branch)
        expr_subset <- expr[, order_tip_1[branch_subset]]
        as.data.frame(t(expr_subset)) %>% mutate(Cell = n():1)
    })
heatmaps[[2]]$Cell <- heatmaps[[2]]$Cell + max(heatmaps[[1]]$Cell) + 2L
heatmaps[[3]]$Cell <- heatmaps[[3]]$Cell + max(heatmaps[[2]]$Cell) + 2L
heatmaps %>%
    bind_rows(.id = 'Branch') %>%
    gather(Gene, Expr, -Branch, -Cell) %>%
    ggplot(aes(Cell, Gene, fill = Expr)) + geom_tile() + scale_fill_viridis_c() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0))
    

#################
# 	dpt DEG analysis
#################
#CD4
library(PseudotimeDE)
#options(mc.cores = 8)
#BPPARAM <- BiocParallel::bpparam()
#BPPARAM$workers <- 8
samplestotake <- 100
library(foreach)

CD4DPT_table = data.frame(cell=colnames(scDataMerged.CD4),pseudotime=dpt$dpt)

#generate dpt estimate for 100
sampleingIndexer <- function(x, LPS_sce) {
  suppressPackageStartupMessages(library(Seurat))
  sample(x = c(1:dim(LPS_sce)[2]), size = 0.8*dim(LPS_sce)[2], replace = FALSE)
}
sampleingFunction <- function(x, sce, DPT_table) {
  suppressPackageStartupMessages(library(SingleCellExperiment))
  suppressPackageStartupMessages(library(slingshot))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(destiny))
  print("packages loaded")

  sce <- sce[, x]
  
	#expr = as.matrix(sce@assays$integrated@data) 
	exprbs <- as.matrix(GetAssayData(sce,assay="integrated",slot="data"))
	
	dmbs <- DiffusionMap(t(exprbs),n_pcs = 50,verbose=TRUE)
	dptbs <- DPT(dmbs)
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(dptbs$dpt))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, as_tibble(DPT_table), by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar2 <- function () 
{
    if (exists(".revoDoParCluster", where = doParallel:::.options)) {
        if (!is.null(doParallel:::.options$.revoDoParCluster)) 
            stopCluster(doParallel:::.options$.revoDoParCluster)
        remove(".revoDoParCluster", envir = doParallel:::.options)
    }
}

#generate samples of 80%
#CD4_subsampled_index <- BiocParallel::bplapply(seq_len(samplestotake) ,, LPS_sce = scDataMerged.CD4, BPPARAM = BPPARAM)
#cl <- parallel::makeCluster(8)
#doParallel::registerDoParallel(cl)
#CD4_subsampled_index = foreach(i = 1:samplestotake, .combine = 'c') %dopar% {
#  sampleingIndexer(i, LPS_sce = scDataMerged.CD4)
#}
#unregister_dopar()
#parallel::stopCluster(cl)
CD4_subsampled_index = list()
for(i in 1:samplestotake){
  CD4_subsampled_index[[i]] = sampleingIndexer(i, LPS_sce = scDataMerged.CD4)
}

#CD4_subsampled_results <- BiocParallel::mclapply(CD4_subsampled_index, , sce = scDataMerged.CD4, DPT_table = CD4DPT_table ,BPPARAM = BPPARAM)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
CD4_subsampled_results = foreach(i = 1:length(CD4_subsampled_index)) %dopar% {
  sampleingFunction(CD4_subsampled_index[[i]], sce = scDataMerged.CD4, DPT_table = CD4DPT_table)
}
parallel::stopCluster(cl)

#linear version
if (FALSE){
	CD4_subsampled_results=list()
	#foreach(i = 1:length(CD4_subsampled_index), .combine = 'c') %do% {
	for (i in 1:length(CD4_subsampled_index)){
		print(paste("processing",i))
	  CD4_subsampled_results[[length(CD4_subsampled_results) + 1]] = sampleingFunction(CD4_subsampled_index[[i]], sce = scDataMerged.CD4, DPT_table = CD4DPT_table)
	}
}

#res <- PseudotimeDE::runPseudotimeDE(gene.vec = tail(head(scDataMerged.CD4@assays$RNA@var.features,n=150),n=50),
expr = as.matrix(scDataMerged.CD4@assays$integrated@data)
res <- PseudotimeDE::runPseudotimeDE(gene.vec = scDataMerged.CD4@assays$RNA@var.features,
                                     ori.tbl = CD4DPT_table,
                                     sub.tbl = CD4_subsampled_results,
                                     mat = expr, ## You can also use a matrix or SeuratObj as the input
                                     model = "nb",
                                     mc.cores = 8,mc.preschedule = TRUE)
                                     #,mc.preschedule =FALSE)
                                     #sub.tbl = LPS_sub_tbl[1:100], ## To save time, use 100 subsamples
                                     
PseudotimeDE::plotCurve(gene.vec = res$gene,
                                        ori.tbl = CD4DPT_table,
                                        mat = scDataMerged.CD4,
                                        model.fit = res$gam.fit)
                                        
#################
# 	dpt DEG analysis - tradeSeq
#################
library(tradeSeq)
library(tidyverse)
counts = scDataMerged.CD4@assays$RNA@counts
cellWeights = data.frame(row.names=colnames(scDataMerged.CD4),pseudotime=rep(1,dim(scDataMerged.CD4)[2]))
pseudotime = data.frame(row.names=colnames(scDataMerged.CD4),pseudotime=dpt$dpt)

#plot dpt with branch
plot(dpt, root=3, 2:3,paths_to=c(1:3), col_by = 'num_cells')
plot(dpt, root=3, 2:3, col_by = 'dpr')
pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/pseudotime_destiny/destiny_CD4_branches.pdf"),width=13,height=8)
plot(dpt, root=2L,      col_by = 'branch')
dev.off()
plot(dpt, root=2L,      col_by = 'dpt')



set.seed(5)
#p = register(bpstart(SnowParam(2)))#bpstart(MulticoreParam(4))
#icMat <- evaluateK(counts = counts,  k = 4:10, 
#                   nGenes = 200, verbose = T,pseudotime = pseudotime, 
#                   cellWeights = cellWeights,parallel =TRUE,BPPARAM=p)
#bpstop(p)
icMat <- evaluateK(counts = counts,  k = 4:10, 
                   nGenes = 200, verbose = T,pseudotime = pseudotime, 
                   cellWeights = cellWeights,parallel =FALSE)
                   
pparm=SnowParam(6,type="SOCK")
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 7, verbose = TRUE, parallel = TRUE, BPPARAM=pparm)
saveRDS(sce,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/GAMfit_tradeseq.rds"))

#no parallel
#sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
#                 nknots = 7, verbose = TRUE)
                 
assoRes <- associationTest(sce)
head(assoRes[order(assoRes$pvalue),])
assogenes = head(rownames(assoRes)[order(assoRes$pvalue)],n=100)
assogenesoverlap = assogenes[which(assogenes %in% names(dpt))]

write.table(assoRes,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/GAMfit_DEGs.tsv"),sep="\t")
assoResHVG = assoRes[which(rownames(assoRes) %in% names(dpt)),]
assoResHVG$padj =  p.adjust(p=assoResHVG$pvalue,method="BH")
write.table(assoResHVG,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/GAMfit_DEGs_HVGs.tsv"),sep="\t")

plotSmoothers(sce, counts, gene = "TPT1")
plotSmoothers(sce, counts, gene = "TIGIT")
plotSmoothers(sce, counts, gene = "FOXP3")
plotSmoothers(sce, counts, gene = "IL32")
plotSmoothers(sce, counts, gene = "FOXP3")
plotGeneCount(curve=sce,  gene = "PARK7")

pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/GAMfit_DEGs_top_50_smoothplot.pdf"),width=3,height=3)
#png(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/GAMfit_DEGs_top_50_smoothplot.png"),width=10,height=15)
#par(mfrow=c(5,2))
assoResHVGplt = assoResHVG[which(!grepl("^RP",rownames(assoResHVG))),]
assoResHVGplt = assoResHVGplt[order(assoResHVGplt$waldStat,decreasing=TRUE),]
for (i in 1:10){
	plot(plotSmoothers(sce, counts, gene = rownames(assoResHVGplt)[i],xlab = paste(rownames(assoResHVGplt)[i],"Pseudotime")))
}
dev.off()
#
plot(dpt, col_by = 'TIGIT', pal = viridis::magma)

heatmaps <-
    na.exclude(unique(dpt$Branch)) %>%
    set_names() %>%
    map(function(branch) {
        branch_subset <- dpt$Branch == branch & !is.na(dpt$Branch)
        exprs_subset <- exprs[, order_tip_1[branch_subset]]
        as.data.frame(t(exprs_subset)) %>% mutate(Cell = n():1)
    })
heatmaps[[2]]$Cell <- heatmaps[[2]]$Cell + max(heatmaps[[1]]$Cell) + 2L
heatmaps[[3]]$Cell <- heatmaps[[3]]$Cell + max(heatmaps[[2]]$Cell) + 2L
heatmaps %>%
    bind_rows(.id = 'Branch') %>%
    gather(Gene, Expr, -Branch, -Cell) %>%
    ggplot(aes(Cell, Gene, fill = Expr)) + geom_tile() + scale_fill_viridis_c() +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0))
    
#CD8
library(tradeSeq)
library(tidyverse)
counts = scDataMerged.CD8@assays$RNA@counts
cellWeights = data.frame(row.names=colnames(scDataMerged.CD8),pseudotime=rep(1,dim(scDataMerged.CD8)[2]))
pseudotime = data.frame(row.names=colnames(scDataMerged.CD8),pseudotime=dptCD8$dpt)

plot(dptCD8)

#################
# 	dpt DEG analysis - tradeSeq + diff curves
#################
plot(dpt, root=2L,      col_by = 'dpt')


library(tradeSeq)
library(tidyverse)
counts = scDataMerged.CD4@assays$RNA@counts
cellWeights = data.frame(row.names=colnames(scDataMerged.CD4),curveTFH=as.character(dpt$Branch) == "3",curveTREG=as.character(dpt$Branch) == "2")
cellWeights$curveTFH[is.na(cellWeights$curveTFH)] = 0
cellWeights$curveTREG[is.na(cellWeights$curveTREG)] = 0
cellWeights$curveTFH[which(as.character(dpt$Branch) == "1")] = 0.5
cellWeights$curveTREG[which(as.character(dpt$Branch) == "1")] = 0.5
pseudotime = data.frame(row.names=colnames(scDataMerged.CD4),curveTFH=dpt$dpt,curveTREG=dpt$dpt)
#remove uninformative cells
keeplist = which(cellWeights$curveTFH != 0 | cellWeights$curveTREG != 0)
counts = counts[,keeplist]
cellWeights = cellWeights[keeplist,]
pseudotime = pseudotime[keeplist,]

#plot dpt with branch
set.seed(5)
#p = register(bpstart(SnowParam(2)))#bpstart(MulticoreParam(4))
#icMat <- evaluateK(counts = counts,  k = 4:10, 
#                   nGenes = 200, verbose = T,pseudotime = pseudotime, 
#                   cellWeights = cellWeights,parallel =TRUE,BPPARAM=p)
#bpstop(p)
pparm=SnowParam(6,type="SOCK")
icMat <- evaluateK(counts = counts,  k = 4:10, 
                   nGenes = 200, verbose = T,pseudotime = pseudotime, 
                   cellWeights = cellWeights,parallel = TRUE, BPPARAM=pparm)
                   
pparm=SnowParam(6,type="SOCK")
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 8, verbose = TRUE, parallel = TRUE, BPPARAM=pparm)
saveRDS(sce,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/GAMfit_tradeseq_withCurves.rds"))


#On subsetted lineages
#TFH
lineageID = "TFH"
LINEAGES_TFH = c("2","3")#,NA)

#fix trajactory
labelcells = rep("black",dim(dpt)[1])
labelcells[which(is.na(dpt$Branch))] = "red"
labelcells[which(dpt$Branch %in% LINEAGES_TFH)] = "red"
labelcells[which(dpt$DC2 < -0.01)] = "black"
#plot(dpt$DC1,dpt$DC2)
plot(dpt, root=2L,      col = labelcells)
labelcellsTF = rep(F,dim(dpt)[1])
labelcellsTF[which(is.na(dpt$Branch))] = T
labelcellsTF[which(dpt$Branch %in% LINEAGES_TFH)] = T
labelcellsTF[which(dpt$DC2 < -0.01)] = F

library(tradeSeq)
library(tidyverse)
counts = scDataMerged.CD4@assays$RNA@counts
cellWeights = data.frame(row.names=colnames(scDataMerged.CD4),curvequery=labelcellsTF)#as.character(dpt$Branch) %in% LINEAGES_TFH)
pseudotime = data.frame(row.names=colnames(scDataMerged.CD4),curvequery=rank(dpt$dpt))
keeplist = which(cellWeights$curvequery == TRUE)
#branch fixing
counts = counts[,keeplist]
cellWeights = cellWeights[keeplist,]
pseudotime = pseudotime[keeplist,]
pparm=SnowParam(6,type="SOCK")
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 7, verbose = TRUE, parallel = TRUE, BPPARAM=pparm)
saveRDS(sce,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_tradeseq_TFH.rds"))
#sce=readRDS(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_tradeseq_TFH.rds"))



#standard func
processDPT <- function(lineageID,LINEAGES_TREGS,dpt,scDataMerged){
	dir.create(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID))
	plot(dpt, root=2L,      col_by = 'branch')
	library(tradeSeq)
	library(tidyverse)
	counts = scDataMerged@assays$RNA@counts
	cellWeights = data.frame(row.names=colnames(scDataMerged),curvequery=as.character(dpt$Branch) %in% LINEAGES_TREGS)
	pseudotime = data.frame(row.names=colnames(scDataMerged),curvequery=dpt$dpt)
	keeplist = which(cellWeights$curvequery == TRUE)
	counts = counts[,keeplist]
	cellWeights = cellWeights[keeplist,]
	pseudotime = pseudotime[keeplist,]
	pparm=SnowParam(6,type="SOCK")
	sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
	                 nknots = 7, verbose = TRUE, parallel = TRUE, BPPARAM=pparm)
	saveRDS(sce,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_tradeseq.rds"))
	#sce=readRDS(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_tradeseq.rds"))
	return (sce)
}
processDPTbySeuratCluster <- function(lineageID,LINEAGES_TREGS,dpt,scDataMerged){
	dir.create(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID))
	#plot(dpt, root=2L,      col_by = 'seurat_clusters')
	library(tradeSeq)
	library(tidyverse)
	counts = scDataMerged@assays$RNA@counts
	cellWeights = data.frame(row.names=colnames(scDataMerged),curvequery=as.character(dpt$seurat_clusters) %in% LINEAGES_TREGS)
	pseudotime = data.frame(row.names=colnames(scDataMerged),curvequery=dpt$dpt)
	keeplist = which(cellWeights$curvequery == TRUE)
	counts = counts[,keeplist]
	cellWeights = cellWeights[keeplist,]
	pseudotime = pseudotime[keeplist,]
	pparm=SnowParam(6,type="SOCK")
	sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
	                 nknots = 7, verbose = TRUE, parallel = TRUE, BPPARAM=pparm)
	saveRDS(sce,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_tradeseq.rds"))
	#sce=readRDS(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_tradeseq.rds"))
	return (sce)
}
subsetCounts <- function(lineageID,LINEAGES_TREGS,dpt,scDataMerged){
	library(tradeSeq)
	library(tidyverse)
	counts = scDataMerged@assays$RNA@counts
	cellWeights = data.frame(row.names=colnames(scDataMerged),curvequery=as.character(dpt$Branch) %in% LINEAGES_TREGS)
	keeplist = which(cellWeights$curvequery == TRUE)
	counts = counts[,keeplist]
	return(counts)
}
subsetCountsbySeuratCluster <- function(lineageID,LINEAGES_TREGS,dpt,scDataMerged){
	library(tradeSeq)
	library(tidyverse)
	counts = scDataMerged@assays$RNA@counts
	cellWeights = data.frame(row.names=colnames(scDataMerged),curvequery=as.character(dpt$seurat_clusters) %in% LINEAGES_TREGS)
	keeplist = which(cellWeights$curvequery == TRUE)
	counts = counts[,keeplist]
	return(counts)
}

#TFH
lineageID = "TFH"
LINEAGES_TREGS = c("2","3")
sce = processDPT(lineageID,LINEAGES_TREGS,dpt,scDataMerged.CD4)
counts = subsetCounts(lineageID,LINEAGES_TREGS,dpt,scDataMerged.CD4)

#TREGS
lineageID = "TREGS"
LINEAGES_TREGS = c("2","1")
sce = processDPT(lineageID,LINEAGES_TREGS,dpt,scDataMerged.CD4)
counts = subsetCounts(lineageID,LINEAGES_TREGS,dpt,scDataMerged.CD4)

#CD8 TEM1
lineageID = "CD8_TEM1"
CLUSTERANA = c("3","1")

sce = processDPT(lineageID,CLUSTERANA,dptCD8,scDataMerged.CD8)
counts = subsetCounts(lineageID,CLUSTERANA,dptCD8,scDataMerged.CD8)

#CD8 TEX
lineageID = "CD8_TEX"
CLUSTERANA = c("2","1")

sce = processDPT(lineageID,CLUSTERANA,dptCD8,scDataMerged.CD8)
counts = subsetCounts(lineageID,CLUSTERANA,dptCD8,scDataMerged.CD8)

#CD8 TEM1 by seurat
lineageID = "CD8_TEM1"
SEURATCLUSTANA = c("7","2")

sce = processDPTbySeuratCluster(lineageID,SEURATCLUSTANA,dptCD8,scDataMerged.CD8)
counts = subsetCountsbySeuratCluster(lineageID,SEURATCLUSTANA,dptCD8,scDataMerged.CD8)

#CD8 TEX by seurat
lineageID = "CD8_TEX"
SEURATCLUSTANA = c("7","8")

sce = processDPTbySeuratCluster(lineageID,SEURATCLUSTANA,dptCD8,scDataMerged.CD8)
counts = subsetCountsbySeuratCluster(lineageID,SEURATCLUSTANA,dptCD8,scDataMerged.CD8)


#################
# 	dpt DEG analysis - slingshot + tradeSeq
#################
dir.create(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID))

if (exists("dptCD8")){ dpt = dptCD8 }
assoRes <- associationTest(sce)
head(assoRes[order(assoRes$pvalue),])
assogenes = head(rownames(assoRes)[order(assoRes$pvalue)],n=100)
assogenesoverlap = assogenes[which(assogenes %in% names(dpt))]

write.table(assoRes,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_DEGs.tsv"),sep="\t")
assoResHVG = assoRes[which(rownames(assoRes) %in% names(dpt)),]
assoResHVG$padj =  p.adjust(p=assoResHVG$pvalue,method="BH")
write.table(assoResHVG,paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_DEGs_HVGs.tsv"),sep="\t")


#plot smoothers plot for top 30
pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/GAMfit_DEGs_top_50_smoothplot.pdf"),width=3,height=3)
#png(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/GAMfit_DEGs_top_50_smoothplot.png"),width=10,height=15)
#par(mfrow=c(5,2))
assoResHVGplt = assoResHVG[which(!grepl("^RP",rownames(assoResHVG))),]
assoResHVGplt = assoResHVGplt[order(assoResHVGplt$waldStat,decreasing=TRUE),]
for (i in 1:30){
	print(i)
	plot(plotSmoothers(sce, counts, gene = rownames(assoResHVGplt)[i],xlab = paste(rownames(assoResHVGplt)[i],"Pseudotime")))
}
dev.off()

#heatmap plot
library(pheatmap)
#mockGenes <-  rownames(assoResHVG)[
#  which(assoResHVG$padj <= 0.05)
#]
mockGenes = rownames(assoResHVGplt)[1:200]
yhatSmooth <- predictSmooth(sce, gene = mockGenes, nPoints = 50, tidy = FALSE)
pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/pheatmap_genes.pdf"),width=10,height=9)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)
#heatSmooth
dev.off()
pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/pheatmap_genes_with_labels.pdf"),width=10,height=30)
pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = TRUE, 
                       show_colnames = FALSE)
dev.off()

#write.table(rownames(yhatSmooth),paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/pheatmap_genes.tsv"),sep="\t")
write.table( heatSmooth$tree_row$labels[heatSmooth$tree_row$order],paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/pheatmap_genes_ordered.tsv"),sep="\t")


#plot the diffusion time 
pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/dpt_by_seurat_clusters.pdf"),width=10,height=9)
plot(dpt,   col_by = 'seurat_clusters')
dev.off()


########################
# k - means clustering
col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)


pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/cheatmap_genes_with_4k.pdf"),width=10,height=9)

HM <-	Heatmap(t(scale(t(yhatSmooth[, 1:50]))), km = 4,
	show_row_names = FALSE, col=col, cluster_columns = FALSE,
	column_title = "Pseudotime clusters",show_column_names =FALSE)
	HM
dev.off()

#HM
pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/cheatmap_genes_with_4k_with_labels.pdf"),width=10,height=40)

HM <-	Heatmap(t(scale(t(yhatSmooth[, 1:50]))), km = 4,
	show_row_names = TRUE, col=col, cluster_columns = FALSE,
	column_title = "Pseudotime clusters",show_column_names =FALSE, row_gap = unit(2, "mm"))
	HM
dev.off()



#k means outside complex heatmap
kclus <- kmeans(t(scale(t(yhatSmooth[, 1:50]))), 4)
col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
head(kclus$cluster, 50)
# Use the split parameter of Heatmap() in order to split the heatmap based on the k-means result:
split <- paste0("Cluster\n", kclus$cluster)
HM <-	Heatmap(t(scale(t(yhatSmooth[, 1:50]))),
	show_row_names = FALSE, col=col, cluster_columns = FALSE,split=split, cluster_row_slices = FALSE,
	column_title = "Pseudotime clusters")
	HM
	
if (lineageID == "TREGS"){
	split = gsub("Cluster\n1","Cluster\n1OLD",split)
	split = gsub("Cluster\n2","Cluster\n1",split)
	split = gsub("Cluster\n1OLD","Cluster\n2",split)
	
	#split = gsub("Cluster\n2","Cluster\n2OLD",split)
	#split = gsub("Cluster\n3","Cluster\n2",split)
	#split = gsub("Cluster\n2OLD","Cluster\n3",split)
}

pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/cheatmap_genes_with_4kmanual.pdf"),width=10,height=9)
HM <-	Heatmap(t(scale(t(yhatSmooth[, 1:50]))),
	show_row_names = FALSE, col=col, cluster_columns = FALSE,split=split, cluster_row_slices = FALSE,
	column_title = "Pseudotime clusters")
	HM
dev.off()

pdf(paste(sep="",projectdir,"/analysis/markers/",mergeID,"/",lineageID,"/cheatmap_genes_with_4kmanual_with_labels.pdf"),width=10,height=40)
HM <-	Heatmap(t(scale(t(yhatSmooth[, 1:50]))),
	show_row_names = TRUE, col=col, cluster_columns = FALSE,split=split, cluster_row_slices = FALSE,
	column_title = "Pseudotime clusters", row_gap = unit(2, "mm"))
	HM
dev.off()