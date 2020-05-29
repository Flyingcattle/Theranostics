
################################################################################
############## Workflow for the C57 Genital Ridges #############################
################################################################################

memory.limit(size = 56000)

library(Seurat)
library(Matrix)
library(dplyr)


# Read in all four input expression matrices

NB <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_skin/cellranger2.1.0_output/0dpp/outs/filtered_gene_bc_matrices/mm10")
E13 <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_skin/cellranger2.1.0_output/13dpc/outs/filtered_gene_bc_matrices/mm10")
E16 <- Read10X(data.dir = "E:/Lab of Germ Cell Bio/10xGenomics_skin/cellranger2.1.0_output/16dpc/outs/filtered_gene_bc_matrices/mm10")


colnames(x = NB) <- paste('0dpp', colnames(x = NB), sep = '_')
colnames(x = E13) <- paste('E13.5', colnames(x = E13), sep = '_')
colnames(x = E16) <- paste('E16.5', colnames(x = E16), sep = '_')


##creat Seurat Object
NB_data <- CreateSeuratObject(raw.data = NB, min.cells = 3, min.genes = 200, project = "odpp")

NB_data <- FilterCells(NB_data, subset.names = "nGene", low.thresholds = 1750)
NB_data <- NormalizeData(NB_data)
NB_data <- FindVariableGenes(NB_data, do.plot = F, display.progress = F)
NB_data <- ScaleData(NB_data)
NB_data@meta.data$tech <- "0dpp"

NB_data <- CreateSeuratObject(raw.data = E13, min.cells = 3, min.genes = 200, project = "E13.5")

E13_data <- FilterCells(E13_data, subset.names = "nGene", low.thresholds = 1750)
E13_data <- NormalizeData(E13_data)
E13_data <- FindVariableGenes(E13_data, do.plot = F, display.progress = F)
E13_data <- ScaleData(E13_data)
E13_data@meta.data$tech <- "E13.5"

E16_data <- CreateSeuratObject(raw.data = E16, min.cells = 3, min.genes = 200, project = "E16.5")

E16_data <- FilterCells(E16_data, subset.names = "nGene", low.thresholds = 1750)
E16_data <- NormalizeData(E16_data)
E16_data <- FindVariableGenes(E16_data, do.plot = F, display.progress = F)
E16_data <- ScaleData(E16_data)
E16_data@meta.data$tech <- "E16.5"


# run CCA 
ob.list <- list(NB_data, E13_data, E16_data)
genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

skin.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 15)        ##not sure here


# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = skin.integrated, reduction.use = "cca", group.by = "tech", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = skin.integrated, features.plot = "CC1", group.by = "tech", point.size.use = 0,
              do.return = TRUE)
plot_grid(p1, p2)


PrintDim(object = skin.integrated, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)


# CC Selection
MetageneBicorPlot(skin.integrated, grouping.var = "tech", dims.eval = 1:15)


DimHeatmap(object = skin.integrated, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:12, do.balanced = TRUE)


##Align the CCA subspaces

skin.integrated <- AlignSubspace(skin.integrated, reduction.type = "cca", grouping.var = "tech", 
                                dims.align = 1:10)

#visualize the aligned CCA and perform an integrated analysis 

p1 <- VlnPlot(object = skin.integrated, features.plot = "ACC1", group.by = "tech", point.size.use = 0,
              do.return = TRUE)
p2 <- VlnPlot(object = skin.integrated, features.plot = "ACC2", group.by = "tech", point.size.use = 0,
              do.return = TRUE)
plot_grid(p1, p2)


## Perform an integrated analysis
# t-SNE and Clustering
skin.integrated <- RunTSNE(skin.integrated, reduction.use = "cca.aligned", dims.use = 1:10, 
                          do.fast = T)


skin.integrated <- FindClusters(skin.integrated, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:10)

# Visualization
p1 <- TSNEPlot(skin.integrated, do.return = T, pt.size = 0.5, group.by = "tech")
p2 <- TSNEPlot(skin.integrated, do.label = F, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)  


#### Find cluster markers 
skin.markers <- FindAllMarkers(object = skin.integrated, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)



##VlnPlot to visualize marked genes
VlnPlot(object = skin.integrated, features.plot = c("Cd34", "Sox9")) 


##########  cell type annotation			
			
## Dermal markers    ## 6.5x3 
FeaturePlot(object = skin.integrated, features.plot = c("Col1a1",'Lum'), cols.use = c("grey", "red"), pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")

dermal <- FeaturePlot(object = skin.integrated, features.plot = c("Col1a1",'Lum'), cols.use = c("grey", "red"), pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Epithelial markers    ## 6.5x3 
FeaturePlot(object = skin.integrated, features.plot = c("Krt14", 'Krt15'), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
epithelial <- FeaturePlot(object = skin.integrated, features.plot = c("Krt14", 'Krt15'), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Dermal condensate markers    ## 6.5x3 
FeaturePlot(object = skin.integrated, features.plot = c("Sox2","Sox18"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
DC <- FeaturePlot(object = skin.integrated, features.plot = c("Sox2","Sox18"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Endothelial markers   ### 6.5x3 
FeaturePlot(object = skin.integrated, features.plot = c("Pecam1",'Kdr'), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
endothelial <- FeaturePlot(object = skin.integrated, features.plot = c("Pecam1",'Kdr'), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Pericytes markers ### 6.5x3     # Acta2= a-SMA
FeaturePlot(object = skin.integrated, features.plot = c("Rgs5","Acta2"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
pericytes <- FeaturePlot(object = skin.integrated, features.plot = c("Rgs5","Acta2"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Melanocytes markers 
FeaturePlot(object = skin.integrated, features.plot = c("Plp1","Fabp7"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
melanocytes <- FeaturePlot(object = skin.integrated, features.plot = c("Plp1","Fabp7"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Immue cell markers 
FeaturePlot(object = skin.integrated, features.plot = c("Cd52","Fcer1g"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
immune <- FeaturePlot(object = skin.integrated, features.plot = c("Cd52","Fcer1g"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Neural markers 
FeaturePlot(object = skin.integrated, features.plot = c("Stmn3","Map2"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
neural <- FeaturePlot(object = skin.integrated, features.plot = c("Stmn3","Map2"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)

## Muscle satellite
FeaturePlot(object = skin.integrated, features.plot = c("Pax7","Myod1"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne")
muscle <- FeaturePlot(object = skin.integrated, features.plot = c("Pax7","Myod1"), cols.use = c("grey", "red") ,pt.size = 0.5, no.axes = TRUE, no.legend = FALSE,
            reduction.use = "tsne", do.return=TRUE)
			
			

####Selecting particular cells and subsetting the Seurat object#################################
########################important#####################################
# use SubsetData to select particular clusters 

dermal.subset <- SubsetData(object = skin.integrated, ident.use = c("Dermal_0", "Dermal_1", "Dermal_2","Dermal_4","Dermal_8","Dermal_condensate_7"))
table(dermal.subset@meta.data$res.0.5)   ##    0    1    2    4    7    8 
                                         ##  2432 2370 2339 1344  813  587 
epithelial.subset <- SubsetData(object = skin.integrated, ident.use = c("Epithelial_3", "Epithelial_5"))
table(epithelial.subset@meta.data$res.0.5)
## rerun the pipeline after subsetting. 

