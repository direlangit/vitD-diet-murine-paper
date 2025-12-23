### Introduction
#===============================================================================
# Title: Buccal Subset Analysis
# Purpose: Subset initial Seurat object based on coarse annotations then perform
#          analysis, functional annotation, clustering, pathway analysis
#===============================================================================

### Packages
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



### Subset the different tissues, buccal
#===============================================================================
buccal_epi <- subset(m_sob.list_mms_buccal,
                     m_sob.list_mms_buccal$coarse_ann2 == 'Epithelia')

buccal_epi_sec <- subset(m_sob.list_mms_buccal,
                         m_sob.list_mms_buccal$coarse_ann2 %in% c('Secretory-Epi'))

buccal_endo <- subset(m_sob.list_mms_buccal,
                      m_sob.list_mms_buccal$coarse_ann2 == c('Endothelia'))

buccal_fibro <- subset(m_sob.list_mms_buccal,
                       m_sob.list_mms_buccal$coarse_ann2 == c('Fibroblasts'))

buccal_obl <- subset(m_sob.list_mms_buccal,
                     m_sob.list_mms_buccal$coarse_ann2 %in% c('Osteoblast-like'))

buccal_myeloid <- subset(m_sob.list_mms_buccal,
                         m_sob.list_mms_buccal$coarse_ann2 == c('Myeloid'))

buccal_skm <- subset(m_sob.list_mms_buccal,
                     m_sob.list_mms_buccal$coarse_ann2 == c('Skeletal_muscle'))
#===============================================================================



### Standard workflow, buccal_epi
#===============================================================================
buccal_epi <- NormalizeData(buccal_epi)
buccal_epi <- FindVariableFeatures(buccal_epi)
buccal_epi <- ScaleData(buccal_epi, vars.to.regress = 'percent.mt')
buccal_epi <- RunPCA(buccal_epi)
ElbowPlot(buccal_epi, ndims = 50)
buccal_epi <- FindNeighbors(buccal_epi,
                            reduction='pca',
                            dims=1:10)
buccal_epi <- FindClusters(buccal_epi,
                           resolution=1,
                           algorithm = 4)
buccal_epi <- RunUMAP(buccal_epi,
                      reduction='pca',
                      dims=1:10)

#===============================================================================

### Findallmarkers, buccal_epi
#===============================================================================
buccal_epi <- JoinLayers(buccal_epi)
Idents(buccal_epi) <- 'seurat_clusters'
buccal_epi_marks <- FindAllMarkers(buccal_epi,
                                   only.pos = T,
                                   min.pct = 0.25
)

top20_be <- buccal_epi_marks %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

DotPlot(buccal_epi,
        features=unique(top20_be$gene),
        cluster.idents=T,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=12),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')
#===============================================================================

### Annotating clusters, buccal_epi
#===============================================================================
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 1] <- 'EP.1'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 2] <- 'EP.2'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 3] <- 'EP.1'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 4] <- 'EP.2'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 5] <- 'EP.1'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 6] <- 'EP.3'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 7] <- 'EP.4'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 8] <- 'EP.5'
buccal_epi@meta.data$subset_ann[buccal_epi@meta.data$seurat_clusters == 9] <- 'EP.6'


umap_buccal_epi <- DimPlot2(
  buccal_epi,
  group.by = 'subset_ann',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'pro_yellow',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)
umap_buccal_epi


cd_buccal_epi <- ClusterDistrBar(
  buccal_epi$condition,
  buccal_epi$subset_ann,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'pro_yellow',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
cd_buccal_epi
#===============================================================================

### Dotplot for subset, buccal_epi
#===============================================================================

Idents(buccal_epi) <- 'subset_ann'
buccal_epi$subset_ann <- factor(buccal_epi$subset_ann, levels=c('EP.6',
                                                                'EP.5',
                                                                'EP.4',
                                                                'EP.3',
                                                                'EP.2',
                                                                'EP.1'
))


buccal_epi_genes <- c(
  'Krt5', 'Krt14', 'Krt15', 'Krt17',
  'Ankrd13b', 'Rhbdd1',
  'Eln', 'Adamts14', 'Dact2',
  'Alox12e', 'Cst6', 'Rgs4',
  'Rhov', 'Lipk',
  'Cxcl14',
  'Myh9', 'Col17a1'
)

DotPlot2(buccal_epi,
         features = unique(buccal_epi_genes),
         flip=T,
         scale_percent = T,
         split.by = 'subset_ann',
         group.by = 'subset_ann',
         split.by.method = 'color',
         split.by.colors = cols_cells_be) +
  theme(axis.text.x = element_text(size=16, face='bold', color='black'),
        axis.text.y = element_text(size=14, face='bold', color='black'),
        legend.key.size = unit(0.1, 'cm')) +
  guides(color='none')


cols_cells_be <- c(
  'EP.1'      = "#db6229",
  'EP.2'       = "#a9816b",
  'EP.3'     = "#6e442a",
  'EP.4'         = "#9e5927",
  'EP.5'   = "#e3a77f",
  'EP.6' = "#dd9835"
)

#===============================================================================

# ORA, buccal_epi
#===============================================================================
buccal_epi <- JoinLayers(buccal_epi)

buccal_epi_de <- FindMarkers(buccal_epi,
                             group.by = 'condition',
                             ident.1 = 'Vit D-',
                             ident.2 = 'Vit D+')

uni_genes <- rownames(m_sob.list_mms_buccal)
de_genes <- rownames(
  subset(buccal_epi_de, p_val_adj <= 0.05)
)

buccal_epi_ora <- enrichGO(de_genes,
                           OrgDb = org.Mm.eg.db,
                           universe = uni_genes,
                           keyType = 'SYMBOL',
                           readable = T,
                           ont = 'ALL',
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.1)


clusterProfiler::dotplot(buccal_epi_ora,
                         showCategory = 20)
#===============================================================================




### Standard workflow, buccal_fibro
#===============================================================================
buccal_fibro <- NormalizeData(buccal_fibro)
buccal_fibro <- FindVariableFeatures(buccal_fibro)
buccal_fibro <- ScaleData(buccal_fibro, vars.to.regress = 'percent.mt')
buccal_fibro <- RunPCA(buccal_fibro)
ElbowPlot(buccal_fibro, ndims = 50)
buccal_fibro <- FindNeighbors(buccal_fibro,
                              reduction='pca',
                              dims=1:10)
buccal_fibro <- FindClusters(buccal_fibro,
                             resolution=1,
                             algorithm = 4)
buccal_fibro <- RunUMAP(buccal_fibro,
                        reduction='pca',
                        dims=1:10)
#===============================================================================

### Findallmarkers Epi, buccal_fibro
#===============================================================================
buccal_fibro <- JoinLayers(buccal_fibro)
Idents(buccal_fibro) <- 'seurat_clusters'
buccal_fibro_marks <- FindAllMarkers(buccal_fibro,
                                     only.pos = T,
                                     min.pct = 0.25
)

top20_fi <- buccal_fibro_marks %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

DotPlot(buccal_fibro,
        features=unique(top20_fi$gene),
        cluster.idents=T,
        scale = F,
        group.by='coarse_ann') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=16),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')

#===============================================================================

### Annotating clusters, buccal_fibro
#===============================================================================

buccal_fibro@meta.data$subset_ann[buccal_fibro@meta.data$seurat_clusters == 1] <- 'FB.1'
buccal_fibro@meta.data$subset_ann[buccal_fibro@meta.data$seurat_clusters == 2] <- 'FB.1'
buccal_fibro@meta.data$subset_ann[buccal_fibro@meta.data$seurat_clusters == 3] <- 'FB.2'
buccal_fibro@meta.data$subset_ann[buccal_fibro@meta.data$seurat_clusters == 4] <- 'FB.3'
buccal_fibro@meta.data$subset_ann[buccal_fibro@meta.data$seurat_clusters == 5] <- 'FB.4'


umap_buccal_fibro <- DimPlot2(
  buccal_fibro,
  group.by = 'subset_ann',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'pro_blue',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)
umap_buccal_fibro

cd_buccal_fibro <- ClusterDistrBar(
  buccal_fibro$condition,
  buccal_fibro$subset_ann,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'pro_blue',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
cd_buccal_fibro

#===============================================================================

### Dotplot for subset, buccal_fibro
#===============================================================================

Idents(buccal_fibro) <- 'subset_ann'

buccal_fibro_save <- buccal_fibro
buccal_fibro$subset_ann <- factor(buccal_fibro$subset_ann, levels=c('FB.4',
                                                                    'FB.2',
                                                                    'FB.3',
                                                                    'FB.1'
))



buccal_fibro_genes <- c(
  'Vim', 'Lum', 'Col1a1', 'Col6a2',
  'Mplkip', 'Polr3a', 'Iscu',
  'Slc1a5', 'Imp3', 'Meg3',
  'Ccdc170', 'Dhdh', 'Pank2'
)

table(buccal_fibro)

DotPlot2(buccal_fibro,
         features = (buccal_fibro_genes),
         flip=T,
         scale_percent = T,
         group.by = 'subset_ann',
         split.by = 'subset_ann',
         split.by.method = 'color',
         split.by.colors = 'pro_blue') +
  theme(axis.text.x = element_text(size=16, face='bold', color='black'),
        axis.text.y = element_text(size=14, face='bold', color='black'),
        legend.key.size = unit(0.1, 'cm')) +
  guides(color='none')

#===============================================================================

### ORA, buccal_fibro
#===============================================================================
buccal_fibro <- JoinLayers(buccal_fibro)

buccal_fibro_de <- FindMarkers(buccal_fibro,
                               group.by = 'condition',
                               ident.1 = 'Vit D-',
                               ident.2 = 'Vit D+')

uni_genes <- rownames(m_sob.list_mms_buccal)
de_genes_fi <- rownames(
  subset(buccal_fibro_de, p_val_adj <= 0.05)
)

buccal_fi_ora <- enrichGO(de_genes_fi,
                          OrgDb = org.Mm.eg.db,
                          universe = uni_genes,
                          keyType = 'SYMBOL',
                          readable = T,
                          ont = 'ALL',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1)

buccal_fi_ora_simp <- simplify(buccal_fi_ora)

clusterProfiler::dotplot(
  buccal_fi_ora_simp,
  showCategory=20)
#===============================================================================




### Standard workflow, buccal_endo
#===============================================================================
buccal_endo <- NormalizeData(buccal_endo)
buccal_endo <- FindVariableFeatures(buccal_endo)
buccal_endo <- ScaleData(buccal_endo, vars.to.regress = 'percent.mt')
buccal_endo <- RunPCA(buccal_endo)
ElbowPlot(buccal_endo, ndims = 50)
buccal_endo <- FindNeighbors(buccal_endo,
                             reduction='pca',
                             dims=1:5)
buccal_endo <- FindClusters(buccal_endo,
                            resolution=1,
                            algorithm = 4)
buccal_endo <- RunUMAP(buccal_endo,
                       reduction='pca',
                       dims=1:5)


pct <- buccal_endo[["pca"]]@stdev / sum(buccal_endo[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1


umap_buccal_fibro <- DimPlot2(
  buccal_endo,
  group.by = 'seurat_clusters',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 5,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'light',
  theme = list(labs(title=NULL, x='UMAP 1', y='UMAP 2'), 
               NoLegend(), theme_umap_arrows(line_length = 25,
                                             text_size = 16,
                                             line_width = 5,
                                             anchor_x = 7,
                                             anchor_y = 7))
)
umap_buccal_fibro

cd_buccal_endo <- ClusterDistrBar(
  buccal_endo$condition,
  buccal_endo$seurat_clusters,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'pro_blue'
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
cd_buccal_endo

#===============================================================================

### Findallmarkers, buccal_endo
#===============================================================================
buccal_endo <- JoinLayers(buccal_endo)
Idents(buccal_endo) <- 'seurat_clusters'
buccal_endo_marks <- FindAllMarkers(buccal_endo,
                                    only.pos = T,
                                    min.pct = 0.25
)


buccal_endo <- JoinLayers(buccal_endo)
Idents(buccal_endo) <- 'subset_ann'
buccal_endo_marks <- FindAllMarkers(buccal_endo,
                                    only.pos = T,
                                    min.pct = 0.25
)

buccal_endo[['RNA']] <- split(buccal_endo[['RNA']], 
                              f=buccal_endo$orig.ident)


top30_endo <- buccal_endo_marks %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)

write.csv(top30_endo, 'top30_endo_marks.csv')

DotPlot(buccal_endo,
        features=unique(c('Cdh5','Pecam1','Prox1','Pdpn','Lyve-1',
                          top30_endo$gene)),
        cluster.idents=T,
        scale = F,
        group.by='subset_ann') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=16),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')

write.csv(top30_endo,
          'buccal_top30_endo.csv')

#===============================================================================

### Annotating clusters, buccal_endo
#===============================================================================

buccal_endo@meta.data$subset_ann[buccal_endo@meta.data$seurat_clusters == 1] <- 'EN.1'
buccal_endo@meta.data$subset_ann[buccal_endo@meta.data$seurat_clusters == 2] <- 'EN.2'
buccal_endo@meta.data$subset_ann[buccal_endo@meta.data$seurat_clusters == 3] <- 'EN.2'

umap_buccal_endo <- DimPlot2(
  buccal_endo,
  group.by = 'subset_ann',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'pro_red',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)
umap_buccal_endo

cd_buccal_endo <- ClusterDistrBar(
  buccal_endo$condition,
  buccal_endo$subset_ann,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'pro_red',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
cd_buccal_endo
#===============================================================================

### Dotplot for subset, buccal_endo
#===============================================================================

buccal_endo_genes <- c(
  'Bcam', 'Cd34', 'Nrp1', 'Prox1', 'Pdpn', 'Lyve1',
  'Lpl', 'Fabp4', 
  'Plvap', 'Il6st'
)

DotPlot2(buccal_endo,
         features = buccal_endo_genes,
         flip=T,
         scale_percent = T,
         group.by = 'subset_ann',
         split.by = 'subset_ann',
         split.by.method = 'color',
         split.by.colors = 'pro_red'
) +
  theme(axis.text.x = element_text(size=16, face='bold', color='black'),
        axis.text.y = element_text(size=14, face='bold', color='black'),
        legend.key.size = unit(0.1, 'cm')) +
  guides(color='none')

#===============================================================================

### ORA, buccal_endo
#===============================================================================
library(DOSE)
buccal_endo <- JoinLayers(buccal_endo)

buccal_endo_de <- FindMarkers(buccal_endo,
                              group.by = 'condition',
                              ident.1 = 'Vit D-',
                              ident.2 = 'Vit D+')

uni_genes_buccal <- rownames(m_sob.list_mms_buccal)
de_genes_endo <- rownames(
  subset(buccal_endo_de, p_val_adj <= 0.05)
)

buccal_endo_ora <- enrichGO(de_genes_endo,
                            OrgDb = org.Mm.eg.db,
                            universe = uni_genes_buccal,
                            keyType = 'SYMBOL',
                            readable = T,
                            ont = 'ALL',
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)

#===============================================================================

### UMAP and cluster distribution bar graph, buccal_endo
#===============================================================================
DimPlot2(
  buccal_endo,
  group.by = 'seurat_clusters',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'light',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)


ClusterDistrBar(
  buccal_endo$condition,
  buccal_endo$seurat_clusters,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'light',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
#===============================================================================



### Standard workflow, buccal_myeloid
#===============================================================================
buccal_myeloid <- NormalizeData(buccal_myeloid)
buccal_myeloid <- FindVariableFeatures(buccal_myeloid)
buccal_myeloid <- ScaleData(buccal_myeloid, vars.to.regress = 'percent.mt')
buccal_epi <- RunPCA(buccal_myeloid)
ElbowPlot(buccal_myeloid, ndims = 50)
buccal_myeloid <- FindNeighbors(buccal_myeloid,
                                reduction='pca',
                                dims=1:10)
buccal_myeloid <- FindClusters(buccal_myeloid,
                               resolution=1,
                               algorithm = 4)
buccal_myeloid <- RunUMAP(buccal_myeloid,
                          reduction='pca',
                          dims=1:10)
#===============================================================================

### Findallmarkers Epi, buccal_myeloid
#===============================================================================
buccal_myeloid <- JoinLayers(buccal_myeloid)
Idents(buccal_myeloid) <- 'seurat_clusters'
buccal_myeloid_marks <- FindAllMarkers(buccal_myeloid,
                                       only.pos = T,
                                       min.pct = 0.25
)


top10_bm <- buccal_myeloid_marks %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


DotPlot(buccal_myeloid,
        features=unique(top10_bm$gene),
        cluster.idents=T,
        scale = F,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=16),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')
#===============================================================================

### ORA buccal_myeloid
#===============================================================================
library(DOSE)
buccal_myeloid <- JoinLayers(buccal_myeloid)

buccal_myeloid_de <- FindMarkers(buccal_myeloid,
                                 group.by = 'condition',
                                 ident.1 = 'Vit D-',
                                 ident.2 = 'Vit D+')

uni_genes_buccal <- rownames(m_sob.list_mms_buccal)
de_genes_myeloid <- rownames(
  subset(buccal_myeloid_de, p_val_adj <= 0.05)
)

buccal_endo_ora <- enrichGO(de_genes_endo,
                            OrgDb = org.Mm.eg.db,
                            universe = uni_genes_buccal,
                            keyType = 'SYMBOL',
                            readable = T,
                            ont = 'ALL',
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.1)




buccal_fi_ora_simp <- simplify(buccal_fi_ora)

clusterProfiler::dotplot(
  buccal_fi_ora_simp,
  showCategory=20)

#===============================================================================

### UMAP and cluster distribution bar graph, buccal_myeloid
#===============================================================================
DimPlot2(
  buccal_myeloid,
  group.by = 'seurat_clusters',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'light',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)


ClusterDistrBar(
  buccal_myeloid$condition,
  buccal_myeloid$seurat_clusters,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'light',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
#===============================================================================



### Standard workflow, buccal_epi_sec
#===============================================================================
buccal_epi_sec <- NormalizeData(buccal_epi_sec)
buccal_epi_sec <- FindVariableFeatures(buccal_epi_sec)
buccal_epi_sec <- ScaleData(buccal_epi_sec, vars.to.regress = 'percent.mt')
buccal_epi_sec <- RunPCA(buccal_epi_sec)
ElbowPlot(buccal_epi_sec, ndims = 50)
buccal_epi_sec <- FindNeighbors(buccal_epi_sec,
                                reduction='pca',
                                dims=1:10)
buccal_epi_sec <- FindClusters(buccal_epi_sec,
                               resolution=1,
                               algorithm = 4)
buccal_epi_sec <- RunUMAP(buccal_epi_sec,
                          reduction='pca',
                          dims=1:10)
#===============================================================================

### Findallmarkers Epi, buccal_epi_sec
#===============================================================================
buccal_epi_sec <- JoinLayers(buccal_epi_sec)
Idents(buccal_epi_sec) <- 'seurat_clusters'
buccal_epi_sec_marks <- FindAllMarkers(buccal_epi_sec,
                                       only.pos = T,
                                       min.pct = 0.25
)


top10_bes <- buccal_epi_sec_marks %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)


DotPlot(buccal_epi_sec,
        features=unique(top10_bes$gene),
        cluster.idents=T,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=16),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')

#===============================================================================

### ORA, buccal_epi_sec
#===============================================================================
buccal_epi_sec <- JoinLayers(buccal_epi_sec)

buccal_epi_sec_de <- FindMarkers(buccal_epi_sec,
                                 group.by = 'condition',
                                 ident.1 = 'Vit D-',
                                 ident.2 = 'Vit D+')

uni_genes <- rownames(m_sob.list_mms_buccal)
de_genes_buccal_epi_sec <- rownames(
  subset(buccal_epi_sec_de, p_val_adj <= 0.05)
)

buccal_epi_sec_ora <- enrichGO(de_genes_buccal_epi_sec,
                               OrgDb = org.Mm.eg.db,
                               universe = uni_genes,
                               keyType = 'SYMBOL',
                               readable = T,
                               ont = 'ALL',
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05)


clusterProfiler::dotplot(buccal_epi_sec_ora)
#===============================================================================

### UMAP and cluster distribution bar graph, buccal_epi_sec
#===============================================================================
DimPlot2(
  buccal_epi_sec,
  group.by = 'seurat_clusters',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'light',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)


ClusterDistrBar(
  buccal_epi_sec$condition,
  buccal_epi_sec$seurat_clusters,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'light',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
#===============================================================================



### Standard workflow, buccal_obl
#===============================================================================
buccal_obl <- NormalizeData(buccal_obl)
buccal_obl <- FindVariableFeatures(buccal_obl)
buccal_obl <- ScaleData(buccal_obl, vars.to.regress = 'percent.mt')
buccal_obl <- RunPCA(buccal_obl, npcs = 48)
ElbowPlot(buccal_obl, ndims = 48)
buccal_obl <- FindNeighbors(buccal_obl,
                            reduction='pca',
                            dims=1:10)
buccal_obl <- FindClusters(buccal_obl,
                           resolution=0.5,
                           algorithm = 4)
buccal_obl <- RunUMAP(buccal_obl,
                      reduction='pca',
                      dims=1:10)
#===============================================================================

### Findallmarkers Epi, buccal_obl
#===============================================================================
buccal_fibro_osteo <- JoinLayers(buccal_fibro_osteo)
Idents(buccal_fibro_osteo) <- 'seurat_clusters'
buccal_fibro_osteo_marks <- FindAllMarkers(buccal_fibro_osteo,
                                           only.pos = T,
                                           min.pct = 0.25
)

top10_be <- buccal_fibro_osteo_marks %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


DotPlot(buccal_fibro_osteo,
        features=unique(top10_be$gene),
        cluster.idents=T,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=16),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')


#===============================================================================

### ORA buccal_fibro_osteo
#===============================================================================
buccal_obl <- JoinLayers(buccal_obl)

buccal_obl_de <- FindMarkers(buccal_obl,
                             group.by = 'condition',
                             ident.1 = 'Vit D-',
                             ident.2 = 'Vit D+')

table(buccal_obl$condition)

uni_genes <- rownames(m_sob.list_mms_buccal)
de_genes_buccal_obl <- rownames(
  subset(buccal_obl_de, p_val_adj <= 0.05)
)

buccal_obl_ora <- enrichGO(de_genes_buccal_obl,
                           OrgDb = org.Mm.eg.db,
                           universe = uni_genes,
                           keyType = 'SYMBOL',
                           readable = T,
                           ont = 'ALL',
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)



dotplot(buccal_obl_ora)

buccal_obl_ora_simp <- simplify(buccal_obl_ora)

clusterProfiler::dotplot(
  buccal_obl_ora_simp,
  showCategory=10)
#===============================================================================

### UMAP and cluster distribution bar graph, buccal_fibro_osteo
#===============================================================================
DimPlot2(
  buccal_obl,
  group.by = 'seurat_clusters',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'light',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)


ClusterDistrBar(
  buccal_obl$condition,
  buccal_obl$seurat_clusters,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'light',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
#===============================================================================




### Standard workflow, buccal_skm
#===============================================================================
buccal_skm <- NormalizeData(buccal_skm)
buccal_skm <- FindVariableFeatures(buccal_skm)
buccal_skm <- ScaleData(buccal_skm, vars.to.regress = 'percent.mt')
buccal_skm <- RunPCA(buccal_skm)
ElbowPlot(buccal_skm, ndims = 50)
buccal_skm <- FindNeighbors(buccal_skm,
                            reduction='pca',
                            dims=1:10)
buccal_skm <- FindClusters(buccal_skm,
                           resolution=1,
                           algorithm = 4)
buccal_skm <- RunUMAP(buccal_skm,
                      reduction='pca',
                      dims=1:10)
#===============================================================================

### Findallmarkers Epi, buccal_skm
#===============================================================================
buccal_skm <- JoinLayers(buccal_skm)
Idents(buccal_skm) <- 'seurat_clusters'
buccal_skm_marks <- FindAllMarkers(buccal_skm,
                                   only.pos = T,
                                   min.pct = 0.25
)

top10_bskm <- buccal_skm_marks %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DotPlot(buccal_skm,
        features=unique(top10_bskm$gene),
        cluster.idents=T,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=16),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')
#===============================================================================

### ORA buccal_skm
#===============================================================================
buccal_skm <- JoinLayers(buccal_skm)

buccal_skm_de <- FindMarkers(buccal_skm,
                             group.by = 'condition',
                             ident.1 = 'Vit D-',
                             ident.2 = 'Vit D+')

uni_genes <- rownames(m_sob.list_mms_buccal)
de_genes_buccal_skm <- rownames(
  subset(buccal_skm_de, p_val_adj <= 0.05)
)

buccal_skm_ora <- enrichGO(de_genes_buccal_skm,
                           OrgDb = org.Mm.eg.db,
                           universe = uni_genes,
                           keyType = 'SYMBOL',
                           readable = T,
                           ont = 'ALL',
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
dotplot(buccal_skm_ora)
#===============================================================================

### UMAP and cluster distribution bar graph, buccal_skm
#===============================================================================
DimPlot2(
  buccal_skm,
  group.by = 'seurat_clusters',
  reduction='umap',
  split.by = 'condition',
  label=T,
  label.size = 3,
  repel=T,
  pt.size = 1.75,
  box=T,
  label.color = 'black',
  cols = 'light',
  theme = list(NoAxes(),labs(title=NULL), NoLegend())
)


ClusterDistrBar(
  buccal_skm$condition,
  buccal_skm$seurat_clusters,
  percent = T,
  flip=F,
  width=0.9,
  border='black',
  cols = 'light',
) + theme(axis.title.y = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, face='bold', angle=45, hjust=1),
          axis.text.y = element_text(size=12, face='bold'),
          legend.text=element_text(size=12, face='bold'),
          legend.title=element_text(size=13, face='bold')) +
  labs(fill='Cluster')
#===============================================================================