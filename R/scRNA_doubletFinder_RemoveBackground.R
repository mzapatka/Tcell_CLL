# Load the DoubletFinder library for doublet detection
	library(DoubletFinder)

# Calculate doublet scores using the doubletScores function from DoubletFinder.
# Parameters:
# - nExp: Number of expected cells per sample (set to 50 here).
# - nDonor: Number of samples or donors in the dataset (set to 10 here).
# - PCs: Principal components to use for analysis (using the first 18 PCs here).
# - use.mito: Boolean indicating whether to use mitochondrial genes for doublet detection (set to FALSE here).
	doub_scores <- doubletScores(seurat_object, 
	                              nExp=50, 
	                              nDonor=10, 
	                              PCs=1:18, 
	                              use.mito=FALSE)
	
# Store the doublet scores in a new metadata column of the Seurat object.
	seurat_object[["doublet_scores"]] <- doub_scores$Score

# Identify cells predicted to be doublets and remove them from the Seurat object.
	doublets <- which(doub_scores$Predicted.Doublet == TRUE)
	seurat_object <- subset(seurat_object, cells = setdiff(colnames(seurat_object), doublets))

# Load the Seurat library for single-cell analysis
	library(Seurat)

# Normalize the data in the Seurat object.
	seurat_object <- NormalizeData(seurat_object)

# Identify variable features using the 'vst' method and retain top 2000 features.
	seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
	seurat_object <- ScaleData(seurat_object)

# Scale the data for downstream analyses.
	
	# Remove ambient RNA
# Define ambient RNA genes based on ribosomal protein gene names.
	ambient_rna_genes <- rownames(seurat_object)[grep("^Rpl[0-9]*|^Rps[0-9]*", rownames(seurat_object))]
	seurat_object <- RemoveBackground(seurat_object, assay = "RNA", cells = NULL, features = ambient_rna_genes)

# Fetch mitochondrial percentage data from the Seurat object for background removal.
	
	data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
	labels <- Idents(seurat_object)

# Remove ambient RNA using the RemoveBackground function in Seurat.
	#meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
	metaDF = read.table(LabelsFile[1],header=T,sep="\t")
# Extract the normalized data matrix from the Seurat object.
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data")

# Retrieve cell labels from the Seurat object to create a metadata dataframe.
	metaMerged = metaDF[colnames(data.input),]
metaDF = read.table(LabelsFile[1], header=TRUE, sep="\t")
	metaMerged = metaMerged[colnames(data.input),]
	