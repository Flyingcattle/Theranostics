

################################################################################
############################# Monocle pipeline  ################################
################################################################################

memory.limit(size = 56000)

library(Seurat)
library(Matrix)
library(dplyr)
library(monocle)


setwd ("E:/Lab of Germ Cell Bio/10xGenomics_skin/@@@Figures_Seurat2.1.0/Figure3")
dermal.subset<-readRDS ("2020.3.7_dermal.subset_res.0.5.rds")


##############  add cellcycle on the trajectory 

head(dermal.subset@meta.data)
exp.mat <- read.table(file = "E:/Lab of Germ Cell Bio/10xGenomics_skin/@@@Figures_Seurat2.1.0/Figure3/Dermal_Monocle/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", 
                      header = TRUE, as.is = TRUE, row.names = 1)
cc.genes <- readLines(con = "E:/Lab of Germ Cell Bio/10xGenomics_skin/@@@Figures_Seurat2.1.0/Figure3/Dermal_Monocle/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

dermal.subset <- CellCycleScoring(object = dermal.subset, s.genes = s.genes, g2m.genes = g2m.genes, 
                           set.ident = TRUE)

RidgePlot(object = dermal.subset, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67"), 
          nCol = 2)

seur_ob<-dermal.subset
## rename seurat object to seur_ob

seur_ob = UpdateSeuratObject(seur_ob)

head (seur_ob@raw.data)
seur_ob@raw.data[1:5,1:5]
seur_ob@data[1:5,1:5]
head (seur_ob@meta.data)

## create CDS from seurat object data
# expression Matrix - select only cells of interest from seurat_ob@data
exprs <- seur_ob@raw.data
cells <- colnames(seur_ob@data)
exprs_mat <- exprs[,match(cells, colnames(exprs))]
exprs_mat[1:5,1:5]

# check number of cells from both
print("raw data dim:")
print(dim(exprs))
print("filtered data dim:")
print(dim(exprs_mat))

## create phenodata from meta data (@data.info for seurat v1)
pheno.data = seur_ob@meta.data

## create feature data\

genes <- data.frame(gene_short_name = rownames(exprs))
head (genes)
rownames(genes) <- rownames(exprs)

# Make CellDataSet object
pd <- new("AnnotatedDataFrame", data=pheno.data)
fd <- new("AnnotatedDataFrame", data=genes)
HSMM_expr_matrix <- exprs_mat

cds <- newCellDataSet(as(HSMM_expr_matrix,"sparseMatrix"),
                      phenoData=pd,
                      featureData=fd,
                      lowerDetectionLimit=0.5,
                      expressionFamily=negbinomial.size())


HSMM <- estimateSizeFactors(cds)
HSMM <- estimateDispersions(HSMM)

## use seurat determined variable genes for ordering
seurat_var_genes = seur_ob@var.genes
head(seur_ob@var.genes)

HSMM_seur_var = setOrderingFilter(HSMM, seurat_var_genes)

plot_ordering_genes(HSMM_seur_var)


## reduce dimensionality with DDRTree, regression on nUMI
HSMM_seur_var <- reduceDimension(HSMM_seur_var, reduction_method="DDRTree",max_components = 2,residualModelFormulaStr = "~nUMI")

## order cells
HSMM_seur_var <- orderCells(HSMM_seur_var)

## print trajectories  6x5
plot_cell_trajectory(HSMM_seur_var, color_by = "Pseudotime",cell_size = 1,show_branch_points = T)
plot_cell_trajectory(HSMM_seur_var, color_by = "tech",cell_size = 1,show_branch_points = T)
plot_cell_trajectory(HSMM_seur_var, color_by = "ori_cluster",cell_size = 1,show_branch_points = T)
plot_cell_trajectory(HSMM_seur_var, color_by = "Phase",cell_size = 1,show_branch_points = T)

### analysis on HSMM object created in monocle scripts ###

HSMM <- HSMM_seur_var

### accessing pheno and feature data ###

pheno.data <- HSMM@phenoData@data
head(pheno.data)
gene.data<-HSMM@featureData@data

plot_cell_trajectory(HSMM_seur_var, color_by = "tech")
plot_cell_trajectory(HSMM_seur_var, color_by = "res.0.5")
plot_cell_trajectory(HSMM_seur_var, color_by = "State")     

### root_state = 1 here
HSMM <- orderCells(HSMM, root_state = 1)

#### Plot cell trajectory
plot_cell_trajectory(HSMM, cell_size=0.5, color_by="State", cell_name_size=3, show_branch_points=T)


###### branch analysis ######

BEAM_res <- BEAM(HSMM, branch_point=2, cores = 4)       


BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

## subset based on qval cut-off  1e-4
subsetted <- subset(BEAM_res, qval < 1e-4)

## plot branched heatmap 

my_col<-colorRampPalette(c(c("#006400","#2E8B57","#FFFF00")))(65)       ## 6x8
plot_genes_branched_heatmap(HSMM[row.names(subsetted),],
                            branch_point=2,
                            num_clusters =4,
                            hmcols = my_col,
                            cores=3,
                            use_gene_short_name = T,
                            show_rownames = F)


