#### Introduction
#===============================================================================
# Title: Gingival Analysis
# Purpose: To perform analysis, clustering, marker identification, and pathway 
#          analysis.
#===============================================================================

#### Packages
#===============================================================================
library(dplyr)
library(Seurat)
library(ggplot2)
library(DropletUtils)
library(clusterProfiler)
library(glmGamPoi)
library(org.Mm.eg.db)
library(enrichplot)
library(BiocManager)
library(here)
library(SeuratExtend)
library(devtools)
set.seed(123)
#===============================================================================

### Standard workflow
#===============================================================================
m_sob.list_mms_ging <- NormalizeData(m_sob.list_mms_ging)
m_sob.list_mms_ging <- FindVariableFeatures(m_sob.list_mms_ging)
m_sob.list_mms_ging <- ScaleData(m_sob.list_mms_ging)
m_sob.list_mms_ging <- RunPCA(m_sob.list_mms_ging)
ElbowPlot(m_sob.list_mms_ging, ndims = 50)
m_sob.list_mms_ging <- FindNeighbors(m_sob.list_mms_ging,
                                     reduction='pca',
                                     dims=1:8)
m_sob.list_mms_ging <- FindClusters(m_sob.list_mms_ging,
                                    resolution=1,
                                    algorithm = 4)

m_sob.list_mms_ging <- RunUMAP(m_sob.list_mms_ging,
                               reduction='pca',
                               dims=1:8)
#===============================================================================

### Differential gene expression, FindAllMarkers
#===============================================================================
m_sob.list_mms_ging <- JoinLayers(m_sob.list_mms_ging)
m_sob.list_mms_ging_marks <- FindAllMarkers(m_sob.list_mms_ging,
                                            only.pos = T,
                                            min.pct = 0.25
)
#===============================================================================

### Gene lists
#===============================================================================
# List of canonical and marker genes from above
list_coarse <- c(
  'Krt5',
  'Krt14',
  'Krt15',
  'Krt17',
  'Kdr',
  'Flt1',
  'Pecam1',
  'Tie1',
  'Col1a1',
  'Lum',
  'Lyz2',
  'Itgam',
  'Acta1',
  'Mb',
  'Tnnc2',
  'Acta2',
  'Myl9'
)
#===============================================================================

### Annotate clusters
#===============================================================================
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 1] <- 'Epithelia'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 2] <- 'Epithelia'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 3] <- 'Epithelia'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 4] <- 'Endothelia'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 5] <- 'Endothelia'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 6] <- 'Endothelia'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 7] <- 'Endothelia'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 8] <- 'Myeloid'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 9] <- 'Fibroblasts'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 10]<- 'Skeletal_muscle'
m_sob.list_mms_ging@meta.data$coarse_ann2[m_sob.list_mms_ging@meta.data$seurat_clusters == 11]<- 'Smooth_muscle'
#===============================================================================

### DEG by condition
#===============================================================================
ging_de <- FindMarkers(m_sob.list_mms_ging,
                       group.by = 'condition',
                       ident.1 = 'Vit D-',
                       ident.2 = 'Vit D+')
#===============================================================================

### Geneset Enrichment Analysis
#===============================================================================
genes_ging <- ging_de$avg_log2FC
names(genes_ging) <- rownames(ging_de)
genes_ging <- na.omit(genes_buccal)
genes_ging <- sort(genes_ging, decreasing=T)

gsea_ging <- gseGO(geneList=genes_ging, 
                           ont ="ALL", 
                           keyType = "SYMBOL", 
                           minGSSize = 3, 
                           maxGSSize = 800, 
                           pvalueCutoff = 0.05,
                           verbose = TRUE,
                           OrgDb = org.Mm.eg.db, 
                           pAdjustMethod = "BH",
                           seed=T,
                           nPermSimple = 10000
)

# Simplify the GSEA terms to reduce redundancy
gse_ging_simp <- simplify(gsea_ging,
                          cutoff = 0.5,
                          by = 'p.adjust',
                          select_fun = min,
                          measure = 'Wang'
)

dot_gse_ging_simp <- enrichplot::dotplot(
  gse_ging_simp, 
  showCategory=10, 
  split='.sign', 
  x='NES',
  font.size=12, 
  label_format=70, 
  title=NULL) +
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=40)
  )
dot_gse_ging_simp + theme(legend.key.size = unit(0.4, 'cm'),
                          legend.text = element_text(size=10))

#===============================================================================

### UMAP generation with color matched cluster distribution bar graph
#===============================================================================

# Setting order
Idents(m_sob.list_mms_ging) <- 'coarse_ann2'
m_sob.list_mms_ging$coarse_ann <- factor(m_sob.list_mms_ging$coarse_ann,
                                   levels=c('Epithelia',
                                            'Endothelia',
                                            'Smooth_muscle',
                                            'Fibroblasts',
                                            'Myeloid',
                                            'Skeletal_muscle'))

# Setting colors
cols_cells_ging <- c(
  'Endothelia'      = "#9e4f53",
  'Epithelia'       = "#ad8251",
  'Fibroblasts'     = "#515b26",
  'Myeloid'         = "#659e77",
  'Smooth_muscle'   = "#718dba",
  'Skeletal_muscle' = "#825189"
)

umap_ging_coarse <- DimPlot2(
  m_sob.list_mms_ging,
  group.by = 'coarse_ann2',
  reduction='umap',
  label=T,
  label.size = 5,
  repel=T,
  pt.size = 1.5,
  box=T,
  label.color = 'black',
  cols=cols_cells_ging,
  theme = list(labs(title=NULL, x='UMAP 1', y='UMAP 2'), 
               NoLegend(), theme_umap_arrows(line_length = 25,
                                             text_size = 16,
                                             line_width = 5,
                                             anchor_x = 7,
                                             anchor_y = 7))
)
umap_ging_coarse

cd_gingiva_coarse <- ClusterDistrBar(
  m_sob.list_mms_ging$condition,
  m_sob.list_mms_ging$coarse_ann2,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols=cols_cells_ging
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cell Type')
cd_gingiva_coarse

#===============================================================================













