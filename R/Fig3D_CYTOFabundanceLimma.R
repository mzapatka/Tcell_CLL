

library(limma)
library(edgeR)
removeControlLN = T


#setwd("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/limmaAnalysis/")
#jj <- read.table("inputDGEallSamples.csv", header = T, sep = ",", row.names = 1)
setwd("CyTOFfcsfileslivecells/CYTOFanalysis_REFER/TcellsSubset/limmaAnalysis/")
jj <- read.table("DafN_Tcells_Table_perSampleClusterCells", header = T, sep = "\t", row.names = 1)
head(jj)
jj1 <- as.matrix(jj)



jj2 <- DGEList(counts = jj1)
condition = as.factor(c(rep("tumorPB", 8), rep("tumorLN", 21), "tumorHodgkinLN",
                        rep("tumorBM",3), rep("tumorAccLN",2), rep("controlLN", 13)))


tissue = as.factor(c(rep("PB", 8), rep("LN",22), rep("BM", 3), rep("LN", 15)))
broad = as.factor(c(rep("T",35), rep("N", 13)))

jj2$samples$condition <- condition
jj2$samples$tissue <- tissue
jj2$samples$broad <- broad


#################### CLL tumor LN VS reactive LN
GT <- jj2$samples$condition
design <- model.matrix(~0+GT)
design


cont.matrix <- makeContrasts(TLNvsRLN=(GTtumorLN-GTcontrolLN),levels=design) #tumorLN-reactiveLN
v <- voom(jj2, design, plot=TRUE)
v
fit <- lmFit(v, design)  
fit2 <- contrasts.fit(fit, cont.matrix)  

fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)

vennDiagram(results)
##############tumorLN-reactiveLN
TLNvsRLN <- topTreat(fit2, coef=1, n=Inf)
head(TLNvsRLN)
write.table(TLNvsRLN, "TLNvsRLN")


####################2 Tumor all VS normal
TN <- jj2$samples$broad
design <- model.matrix(~0+TN)
design


cont.matrix <- makeContrasts(TvsN=(TNT-TNN),levels=design) #tumorAll-reactive
v <- voom(jj2, design, plot=TRUE)
v
fit <- lmFit(v, design)  
fit2 <- contrasts.fit(fit, cont.matrix)  

fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)

vennDiagram(results)
###
TvsN <- topTreat(fit2, coef=1, n=Inf)
head(TvsN)
write.table(TvsN, "TvsN")
#######

####################2 Tumor all+tissue VS normal+tissue=Des3 (tumorPB+tumorBM vs tumorLN)
Des3 <- paste(jj2$samples$broad, jj2$samples$tissue, sep = ".")
design <- model.matrix(~0+Des3)
design


cont.matrix <- makeContrasts(LNvsPBBM=((Des3T.LN+Des3T.BM+Des3T.PB)-Des3N.LN),levels=design) 
v <- voom(jj2, design, plot=TRUE)
v
fit <- lmFit(v, design)  
fit2 <- contrasts.fit(fit, cont.matrix)  

fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)

vennDiagram(results)
##############tumorLN-reactiveLN
LNvsPBBMextra <- topTreat(fit2, coef=1, n=Inf)
head(LNvsPBBM)
write.table(LNvsPBBM, "LNvsPBBM")
#######