library(HDCytoData)
library('tidyverse')
library('readxl')
library('flowCore')
library('Cairo')
library("CATALYST")
library(corrplot)
library(ComplexHeatmap)
library(Seurat)

####  load object for T cells for all the plots below
####   wherever analysis is on all the cells in the dataset, it is mentioned there specifically

#module load R/4.0.3-foss-2020a


daf2 <- readRDS("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/TcellsDaf2.rds")

#############################  UMAP DENSITY PLOT FACET BY CONDITION  ##################################################################


#condition_colors = c("tumorLN"="#0570b0","tumorPB"="#54278f","tumorBM"="#238443", "controlLN" = "#b30000")
condition_colors = c("tumorLN"="#068FDE","tumorPB"="#5F2BA4","tumorBM"="#268446", "controlLN" = "#B30000")
color_sampling_cnt = c(11,9,8,8)

plotDR(daf2, color_by = "condition", facet_by = 'condition') + scale_color_manual(values = condition_colors) + 
  geom_point(size=0.3) + geom_density_2d(size=0.20, colour="black") + ylim(-10, 12) + xlim(-10, 12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

Separateplots = levels(daf2$condition)
library(ggrastr)

for (i in 1:length(Separateplots)){
  #daf2subset = filterSCE(daf2, condition == Separateplots[i])
  #getPalette = colorRampPalette(brewer.pal(9, "Blues"))
  #fillcol = getPalette(11)
  getPalette = colorRampPalette(c("white",condition_colors[Separateplots[i]]))
  #fillcol = getPalette(12)[1:12]
  
  fillcollist=list()
  fillcollist[[1]] = "#FFFFFF00"
  pickColors = color_sampling_cnt[i]
  offset = ceiling(pickColors  / 2)
  fillcolCandidates = getPalette(pickColors + offset)[1:(pickColors+offset)]
  accel = pickColors
  for (j in 2:pickColors){
    accel = round(j/2)#accel - 2
    print(j+accel)
    fillcollist[[j]] = alpha(fillcolCandidates[j+accel],0.5 + (j / pickColors) * 0.5)
  }
  fillcol = do.call(c,fillcollist)
  
  daf2subset = daf2[,which(daf2$condition == Separateplots[i])]
  
  plot<-plotDR(daf2subset, color_by = "condition") + scale_color_manual(values = "gray") + 
    geom_point(size=0.2) + geom_density_2d_filled(size=0.2, colour = "black")  + scale_fill_manual(values  = fillcol) + ylim(-10, 12) + xlim(-10, 12) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
  #plot
  ggplot2::ggsave(filename = paste(sep="","CyToF/density/cytof_plots_",Separateplots[i],".png"),width=6,height=4,dpi=600, plot=plot)
  ggplot2::ggsave(filename = paste(sep="","CyToF/density/cytof_plots_",Separateplots[i],".pdf"),width=6,height=4,dpi=600, plot=plot)
  
  
}
#geom_density_2d_filled(size=0.20, colour="red")
