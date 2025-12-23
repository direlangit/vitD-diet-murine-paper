#### Introduction
#===============================================================================
# Title: Buccal Analysis
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
m_sob.list_mms_buccal <- NormalizeData(m_sob.list_mms_buccal)
m_sob.list_mms_buccal <- FindVariableFeatures(m_sob.list_mms_buccal)
m_sob.list_mms_buccal <- ScaleData(m_sob.list_mms_buccal)
m_sob.list_mms_buccal <- RunPCA(m_sob.list_mms_buccal)
ElbowPlot(m_sob.list_mms_buccal, ndims = 50)
m_sob.list_mms_buccal <- FindNeighbors(m_sob.list_mms_buccal,
                                       reduction='pca',
                                       dims=1:10)
m_sob.list_mms_buccal <- FindClusters(m_sob.list_mms_buccal,
                                      resolution=1,
                                      algorithm = 4)
m_sob.list_mms_buccal <- RunUMAP(m_sob.list_mms_buccal,
                                 reduction='pca',
                                 dims=1:10)
#===============================================================================

### Differential gene expression, FindAllMarkers
#===============================================================================
m_sob.list_mms_buccal <- JoinLayers(m_sob.list_mms_buccal)
m_sob.list_mms_buccal_marks <- FindAllMarkers(
  m_sob.list_mms_buccal,
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
  'Cdh1',
  'Epcam',
  'Col17a1',
  'Dhcr24',
  'Slc27a4',
  'Aco1',
  'Cidea',
  'Lyz2',
  'Itgam',
  'Col1a1',
  'Lum',
  'Acta1',
  'Mb',
  'Tnnc2',
  'Kdr',
  'Flt1',
  'Pecam1',
  'Tie1',
  'Bmp7',
  'Bmp4',
  'Bmp3',
  'Postn',
  'Fgfr1',
  'Mmp11',
  'Col23a1'
)
#===============================================================================

### Annotate clusters
#===============================================================================
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 1] <- 'Epithelia'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 2] <- 'Epithelia'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 3] <- 'Fibroblasts'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 4] <- 'Endothelia'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 5] <- 'Epithelia'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 6] <- 'Secretory-Epi'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 7] <- 'Skeletal_muscle'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 8] <- 'Myeloid'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 9] <- 'Osteoblast-like'
m_sob.list_mms_buccal@meta.data$coarse_ann[m_sob.list_mms_buccal@meta.data$seurat_clusters == 10]<- 'Epithelia'
#===============================================================================

### DEG by condition
#===============================================================================
buccal_de <- FindMarkers(m_sob.list_mms_buccal,
                         group.by = 'condition',
                         ident.1 = 'Vit D-',
                         ident.2 = 'Vit D+')
#===============================================================================

### Gene Set Enrichment Analysis
#===============================================================================
# Create rank list
genes_buccal <- buccal_de$avg_log2FC
names(genes_buccal) <- rownames(buccal_de)
genes_buccal <- na.omit(genes_buccal)
genes_buccal <- sort(genes_buccal, decreasing=T)


gsea_buccal <- gseGO(geneList=genes_buccal, 
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
gse_buccal_simp <- simplify(genes_buccal_nodown,
                            cutoff = 0.5,
                            by = 'p.adjust',
                            select_fun = min,
                            measure = 'Wang'
)

# Dotplot of simplified GSEA result
dot_gse_buccal_simp <- enrichplot::dotplot(
  gse_buccal_simp, 
  showCategory=10, 
  split='.sign', 
  x='NES',
  font.size=12, 
  label_format=70, 
  title=NULL) +
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=40)
  )
dot_gse_buccal_simp + theme(legend.key.size = unit(0.4, 'cm'),
                            legend.text = element_text(size=10))
#===============================================================================

### UMAP generation with color matched cluster distribution bar graph
#===============================================================================

# Setting order
m_sob.list_mms_buccal$coarse_ann <- factor(
  m_sob.list_mms_buccal$coarse_ann, 
  levels=c(
    'Epithelia',
    'Secretory-Epi',
    'Myeloid',
    'Fibroblasts',
    'Skeletal_muscle',
    'Endothelia',
    'Osteoblast-like')
  )


#UMAP
DimPlot2(
  m_sob.list_mms_buccal,
  group.by = 'coarse_ann',
  reduction='umap',
  label=T,
  label.size = 5,
  repel=T,
  pt.size = 1.5,
  box=T,
  label.color = 'black',
  cols='default',
  theme = list(labs(title=NULL, x='UMAP 1', y='UMAP 2'), 
               NoLegend(), theme_umap_arrows(line_length = 25,
                                             text_size = 16,
                                             line_width = 5,
                                             anchor_x = 7,
                                             anchor_y = 7))
)


# CD bar color matched to UMAP
ClusterDistrBar(
  m_sob.list_mms_buccal$condition,
  m_sob.list_mms_buccal$coarse_ann,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols='default'
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cell Type')


#===============================================================================

