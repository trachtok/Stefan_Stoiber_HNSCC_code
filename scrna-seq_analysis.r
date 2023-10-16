#' scRNA-seq analysis of GSE181919 dataset for Stefan Stoiber, 17.5.2023
#' @author Karolina Trachtova
#' 

library(Seurat)
library(ComplexHeatmap)
library(gridExtra)
library(cowplot)
library(gplots)

project_dir <- "~/stefan_scRNA-seq_re-analysis"

################################################################################
# Load input data
################################################################################
# counts downloaded from GEO
# to run the following code, first download files (GSE181919_UMI_counts.txt and GSE181919_Barcode_metadata.txt)
# from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919
# place these files to the ./data folder

# unzip counts matrix -> might take some time or do it manually in the terminal
scrnaseq_counts <- read.table(paste0(project_dir,"/data/GSE181919_UMI_counts.txt"), sep="\t", header=T)

# metadata, downloaded from GEO
metadata <- read.table(paste0(project_dir,"/data/GSE181919_Barcode_metadata.txt"), sep="\t", heade=T)

#  Tumor Initiating Cells -> list of marker genes, derived from the picture from Schramek Lab
TICgenes <- c("CD44", "NGFR", "TP63", "ETS2", "ITGA6", "WNT10A", "MYC", "JUN", "FOS", "KRT5", "KRT14", "KRT10", "KRT1", "MKI67", "CCNA2", "CCNB1", "CCNB2", "KLF4", "PDK1")
SHHgenes <- c("BMP2", "BMP4", "BMP5", "BMP6", "BMP7", "BMP8A", "BMP8B", "BTRC", "CSNK1A1", "CSNK1A1L", "CSNK1D", "CSNK1E", "CSNK1G1", "CSNK1G2", "CSNK1G3", "DHH", "FBXW11", "GAS1", "GLI1", "GLI2", "GLI3", "GSK3B", "HHIP", "IHH", "LRP2", "PRKACA", "PRKACB", "PRKACG", "PRKX", "PTCH1", "PTCH2", "RAB23", "SHH", "SMO" , "STK36" , "SUFU" , "WNT1" , "WNT10A" , "WNT10B" , "WNT11" , "WNT16", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A", "WNT8B", "WNT9A", "WNT9B", "ZIC2")

# check which SHH genes are not in the GSE181919 dataset
SHHgenes[!SHHgenes %in% rownames(scrnaseq)]
# [1] "BMP8A"  "IHH"    "LRP2"   "PRKACG" "SHH"    "WNT1"   "WNT8A"  "WNT8B"  "WNT9B" 

################################################################################
# Data preparation
################################################################################

# filter out tissue specimens we do not want -> keep just primary cancer cells ("CA")
metadata_CA <- metadata[metadata$tissue.type %in% c("CA"),]
rownames(metadata_CA) <- gsub("-","\\.",rownames(metadata_CA))
dim(metadata_CA) # 23088 cells are "CA"

# subset also the initial expression matrix
scrnaseq_CA <- scrnaseq[,colnames(scrnaseq) %in% rownames(metadata_CA)]
dim(scrnaseq_CA)

# create Seurat object
scrnaseq_obj <- CreateSeuratObject(scrnaseq_CA, project = "SeuratProject", assay = "RNA",
                                   min.cells = 0, min.features = 0, names.field = 1,
                                   names.delim = ".", meta.data = metadata_CA)

################################################################################
# FirstQC + pre-filtering
################################################################################

# create slot for MT RNA
scrnaseq_obj[["percent.mt"]] <- PercentageFeatureSet(scrnaseq_obj, pattern = "^MT-")
VlnPlot(scrnaseq_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(scrnaseq_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scrnaseq_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Pre-filtering -> if I wanted to filter out more cells it would be done like this:
#hnscc <- subset(scrnaseq_obj, subset = percent.mt < 10)

################################################################################
# Normalization
################################################################################

hnscc_allFeatures <- SCTransform(scrnaseq_obj, vars.to.regress = "percent.mt", verbose = FALSE, return.only.var.genes = F)

################################################################################
# UMAP + heatmap plotting
################################################################################

# calculate SHH enrichment
hnscc_allFeatures <- AddModuleScore(hnscc_allFeatures,
                                    features = list(SHHgenes),
                                    name="SHH_enriched", ctrl = 10)

dim <- 10
res <- 0.05

hnscc_sb <- subset(x = hnscc_allFeatures, features = TICgenes)
hnscc_sb <- RunPCA(hnscc_sb, verbose = FALSE, approx=FALSE)
    
hnscc_sb <- RunUMAP(hnscc_sb, dims = 1:dim, verbose = FALSE)
hnscc_sb <- FindNeighbors(hnscc_sb, dims = 1:dim, verbose = FALSE)
hnscc_sb <- FindClusters(hnscc_sb, verbose = FALSE, resolution = res)

# Plot UMAP of TIC genes only
p1 <- DimPlot(hnscc_sb, label = TRUE) + ggtitle(paste0("dim 1:", dim, ", resolution=",res))
    
df_mean <- list()
    for(i in 0:(length(levels(hnscc_sb@meta.data$seurat_clusters))-1)){
      clstr <- hnscc_sb@assays$SCT[,colnames(hnscc_sb@assays$SCT) %in% rownames(hnscc_sb@meta.data[hnscc_sb@meta.data$seurat_clusters == i,])]
      tmp <- as.data.frame(rowMeans(clstr))
      colnames(tmp) <- paste0("c",i)
      df_mean <- c(df_mean, list(tmp))
    }

df <- do.call("cbind", df_mean)

# plot heatmap of clusters expression
p2 <- grid.grabExpr(draw(ComplexHeatmap::pheatmap(df, scale="row",
                                                  cluster_cols = F,
                                                  cluster_rows = F,
                                                  color = bluered(201),
                                                  border_color = NA,
                                                  gaps_row = c(13),
                                                  column_names_side = c("top"),  angle_col = c("0"),
                                                  main = "TIC genes (average)")))
    
    
print(plot_grid(p1, p2, nrow=1))

# plot SHH enrichment
p3 <- FeaturePlot(hnscc_sb,
                  features = "SHH_enriched1", label = TRUE, repel = TRUE,
                  order = T,
                  min.cutoff = 0) 
    
print(plot_grid(p1, p2, p3, nrow=1))
