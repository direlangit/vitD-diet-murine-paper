### Introduction
#===============================================================================
# Title: Gingival Subset Analysis
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



### Subset the different tissues, ging
#===============================================================================
ging_epi <- subset(m_sob.list_mms_ging,
                   m_sob.list_mms_ging$coarse_ann == 'Epithelia')

ging_endo <- subset(m_sob.list_mms_ging,
                    m_sob.list_mms_ging$coarse_ann == c('Endothelia'))

ging_fibro <- subset(m_sob.list_mms_ging,
                     m_sob.list_mms_ging$coarse_ann == c('Fibroblasts'))

ging_myeloid <- subset(m_sob.list_mms_ging,
                       m_sob.list_mms_ging$coarse_ann == c('Myeloid'))

ging_skm <- subset(m_sob.list_mms_ging,
                   m_sob.list_mms_ging$coarse_ann == c('Skeletal_muscle'))

ging_sm <- subset(m_sob.list_mms_ging,
                  m_sob.list_mms_ging$coarse_ann == c('Smooth_muscle'))
#===============================================================================



### Standard workflow, ging_epi
#===============================================================================
ging_epi <- NormalizeData(ging_epi)
ging_epi <- FindVariableFeatures(ging_epi)
ging_epi <- ScaleData(ging_epi, vars.to.regress = 'percent.mt')
ging_epi <- RunPCA(ging_epi)
ElbowPlot(ging_epi, ndims = 50)
ging_epi <- FindNeighbors(ging_epi,
                          reduction='pca',
                          dims=1:5)
ging_epi <- FindClusters(ging_epi,
                         resolution=,
                         algorithm = 4)
ging_epi <- RunUMAP(ging_epi,
                    reduction='pca',
                    dims=1:5)
#===============================================================================

### Findallmarkers Epi, ging_epi
#===============================================================================
ging_epi <- JoinLayers(ging_epi)
Idents(ging_epi) <- 'seurat_clusters'
ging_epi_marks <- FindAllMarkers(ging_epi,
                                 only.pos = T,
                                 min.pct = 0.25
)

top10_ge <- ging_epi_marks %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)


DotPlot(ging_epi,
        features=unique(top10_ge$gene),
        cluster.idents=T,
        scale = F,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=12),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')
#===============================================================================

### Annotating clusters, ging_epi
#===============================================================================
ging_epi@meta.data$subset_ann[ging_epi@meta.data$seurat_clusters == 1] <- 'EP.1'
ging_epi@meta.data$subset_ann[ging_epi@meta.data$seurat_clusters == 2] <- 'EP.2'
ging_epi@meta.data$subset_ann[ging_epi@meta.data$seurat_clusters == 3] <- 'EP.2'
ging_epi@meta.data$subset_ann[ging_epi@meta.data$seurat_clusters == 4] <- 'EP.2'
ging_epi@meta.data$subset_ann[ging_epi@meta.data$seurat_clusters == 5] <- 'EP.3'
ging_epi@meta.data$subset_ann[ging_epi@meta.data$seurat_clusters == 6] <- 'EP.4'
ging_epi@meta.data$subset_ann[ging_epi@meta.data$seurat_clusters == 7] <- 'EP.5'

umap_ging_epi <- DimPlot2(
  ging_epi,
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
umap_ging_epi

cd_ging_epi <- ClusterDistrBar(
  ging_epi$condition,
  ging_epi$subset_ann,
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
cd_ging_epi
#===============================================================================

### Dotplot for subset, ging_epi
#===============================================================================
Idents(ging_epi) <- 'subset_ann'
ging_epi$subset_ann <- factor(ging_epi$subset_ann,
                              levels=c(
                                'EP.5',
                                'EP.1',
                                'EP.3',
                                'EP.4',
                                'EP.2'))

ging_epi_list <- c(
  'Krt5',
  'Krt17',
  'Krt13',
  'Krt4',
  'Mt4',
  'Jun',
  'Fos',
  'Egr1',
  'Atf3',
  'Mki67',
  'Top2a',
  'Col16a1',
  'Tnc',
  'Krt75')

DotPlot2(ging_epi,
         features = ging_epi_list,
         flip=T,
         scale_percent = T,
         group.by = 'subset_ann',
         split.by = 'subset_ann',
         split.by.method = 'color',
         split.by.colors = 'pro_yellow') +
  theme(axis.text.x = element_text(size=16, face='bold', color='black'),
        axis.text.y = element_text(size=14, face='bold', color='black'),
        legend.key.size = unit(0.1, 'cm')) +
  guides(color='none')
#===============================================================================

### ORA, ging_epi
#===============================================================================
ging_epi <- JoinLayers(ging_epi)

ging_epi_de <- FindMarkers(ging_epi,
                           group.by = 'condition',
                           ident.1 = 'Vit D-',
                           ident.2 = 'Vit D+')

uni_genes_ging <- rownames(m_sob.list_mms_ging)
de_genes_ging_epi <- rownames(
  subset(ging_epi_de, p_val_adj <= 0.05)
)

ging_epi_ora <- enrichGO(de_genes_ging_epi,
                         OrgDb = org.Mm.eg.db,
                         universe = uni_genes_ging,
                         keyType = 'SYMBOL',
                         readable = T,
                         ont = 'ALL',
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1)

ging_epi_ora_simp <- simplify(ging_epi_ora,
                              cutoff=0.4,
                              by='p.adjust',
                              select_fun = min,
                              measure='Wang')

clusterProfiler::dotplot(ging_epi_ora_simp,
                         showCategory = 10)


View(buccal_epi_ora@result)
#===============================================================================



### Standard workflow, ging_fibro
#===============================================================================
ging_fibro <- NormalizeData(ging_fibro)
ging_fibro <- FindVariableFeatures(ging_fibro)
ging_fibro <- ScaleData(ging_fibro, vars.to.regress = 'percent.mt')
ging_fibro <- RunPCA(ging_fibro)
ElbowPlot(ging_fibro, ndims = 50)
ging_fibro <- FindNeighbors(ging_fibro,
                            reduction='pca',
                            dims=1:10)
ging_fibro <- FindClusters(ging_fibro,
                           resolution=1,
                           algorithm = 4)
ging_fibro <- RunUMAP(ging_fibro,
                      reduction='pca',
                      dims=1:10)
#===============================================================================

### Findallmarkers, ging_fibro
#===============================================================================
ging_fibro <- JoinLayers(ging_fibro)
Idents(ging_fibro) <- 'seurat_clusters'
ging_fibro_marks <- FindAllMarkers(ging_fibro,
                                   only.pos = T,
                                   min.pct = 0.25
)

top20_gf <- ging_fibro_marks %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

DotPlot(ging_fibro,
        features=unique(top20_gf$gene),
        cluster.idents=T,
        scale = F,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=12),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')
#===============================================================================

### Annotating clusters, ging_fibro
#===============================================================================

ging_fibro@meta.data$subset_ann[ging_fibro@meta.data$seurat_clusters == 1] <- 'Fb.1'
ging_fibro@meta.data$subset_ann[ging_fibro@meta.data$seurat_clusters == 2] <- 'Fb.2'


umap_ging_fibro<- DimPlot2(
  ging_fibro,
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
umap_ging_fibro

cd_ging_fibro <- ClusterDistrBar(
  ging_fibro$condition,
  ging_fibro$subset_ann,
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
cd_ging_fibro
#===============================================================================



### Standard workflow, ging_endo
#===============================================================================
ging_endo <- NormalizeData(ging_endo)
ging_endo <- FindVariableFeatures(ging_endo)
ging_endo <- ScaleData(ging_endo, vars.to.regress = 'percent.mt')
ging_endo <- RunPCA(ging_endo)
ElbowPlot(ging_endo, ndims = 50)
ging_endo <- FindNeighbors(ging_endo,
                           reduction='pca',
                           dims=1:5)
ging_endo <- FindClusters(ging_endo,
                          resolution=1,
                          algorithm = 4)
ging_endo <- RunUMAP(ging_endo,
                     reduction='pca',
                     dims=1:5)
#===============================================================================

### Findallmarkers, ging_endo
#===============================================================================
ging_endo <- JoinLayers(ging_endo)
Idents(ging_endo) <- 'seurat_clusters'
ging_endo_marks <- FindAllMarkers(ging_endo,
                                  only.pos = T,
                                  min.pct = 0.25
)

top10_gen <- ging_endo_marks %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


DotPlot(ging_endo,
        features=unique(top10_gen$gene),
        cluster.idents=T,
        scale = F,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=12),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')
#===============================================================================

### Annotating clusters, ging_endo
#===============================================================================
ging_endo@meta.data$subset_ann[ging_endo@meta.data$seurat_clusters == 1] <- 'EN.1'
ging_endo@meta.data$subset_ann[ging_endo@meta.data$seurat_clusters == 2] <- 'EN.2'
ging_endo@meta.data$subset_ann[ging_endo@meta.data$seurat_clusters == 3] <- 'EN.1'
ging_endo@meta.data$subset_ann[ging_endo@meta.data$seurat_clusters == 4] <- 'EN.1'
ging_endo@meta.data$subset_ann[ging_endo@meta.data$seurat_clusters == 5] <- 'EN.4'
ging_endo@meta.data$subset_ann[ging_endo@meta.data$seurat_clusters == 6] <- 'EN.3'

umap_ging_endo<- DimPlot2(
  ging_endo,
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
umap_ging_endo

cd_ging_endo <- ClusterDistrBar(
  ging_endo$condition,
  ging_endo$subset_ann,
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
cd_ging_endo
#===============================================================================

### Dotplot for subset, ging_endo
#===============================================================================
Idents(ging_epi) <- 'subset_ann'

test <- ging_epi
ging_epi_save <- ging_epi

ging_epi$subset_ann <- factor(ging_epi$subset_ann, levels=c('EP.4',
                                                            'EP.1',
                                                            'EP.3',
                                                            'EP.2'
))

ging_endo_genes <- c(
  'Kdr',
  'Git2','Vezf1','Plcb1',
  'Sema3g', 'Col8a1', 'Igf2',
  'Plvap', 'Selp', 'Fgl2'
)


DotPlot2(ging_endo,
         features = ging_endo_genes,
         flip=T,
         scale_percent = T,
         group.by = 'subset_ann',
         split.by = 'subset_ann',
         split.by.method = 'color',
         split.by.colors = 'pro_red') +
  theme(axis.text.x = element_text(size=16, face='bold', color='black'),
        axis.text.y = element_text(size=14, face='bold', color='black'),
        legend.key.size = unit(0.1, 'cm')) +
  guides(color='none')

#===============================================================================

### ORA, ging_endo
#===============================================================================
ging_endo <- JoinLayers(ging_endo)
Idents (ging_endo) <- 'condition'
ging_endo_de <- FindMarkers(ging_endo,
                            group.by = 'condition',
                            ident.1 = 'Vit D-',
                            ident.2 = 'Vit D+')

uni_genes_ging <- rownames(m_sob.list_mms_ging)
de_genes_ging_endo <- rownames(
  subset(ging_endo_de, p_val_adj <= 0.05)
)

ging_endo_ora <- enrichGO(de_genes_ging_endo,
                          OrgDb = org.Mm.eg.db,
                          universe = uni_genes_ging,
                          keyType = 'SYMBOL',
                          readable = T,
                          ont = 'ALL',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.1)

ging_endo_ora_simp <- simplify(ging_endo_ora,
                               cutoff=0.5,
                               by='p.adjust',
                               select_fun = min,
                               measure='Wang')

clusterProfiler::dotplot(ging_endo_ora_simp,
                         showCategory = 10)



View(ging_endo_ora_simp@result)
#===============================================================================



### Standard workflow, ging_skm
#===============================================================================
ging_skm <- NormalizeData(ging_skm)
ging_skm <- FindVariableFeatures(ging_skm)
ging_skm <- ScaleData(ging_skm, vars.to.regress = 'percent.mt')
ging_skm <- RunPCA(ging_skm)
ElbowPlot(ging_skm, ndims = 50)
ging_skm <- FindNeighbors(ging_skm,
                          reduction='pca',
                          dims=1:10)
ging_skm <- FindClusters(ging_skm,
                         resolution=1,
                         algorithm = 4)
ging_skm <- RunUMAP(ging_skm,
                    reduction='pca',
                    dims=1:10)


umap_ging_skm<- DimPlot2(
  ging_skm,
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
umap_ging_skm

cd_ging_skm<- ClusterDistrBar(
  ging_skm$condition,
  ging_skm$seurat_clusters,
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
cd_ging_skm
#===============================================================================

### Findallmarkers Epi, ging_skm
#===============================================================================
ging_endo <- JoinLayers(ging_endo)
Idents(ging_endo) <- 'seurat_clusters'
ging_endo_marks <- FindAllMarkers(ging_endo,
                                  only.pos = T,
                                  min.pct = 0.25
)

top10_gen <- ging_endo_marks %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


DotPlot(ging_endo,
        features=unique(top10_gen$gene),
        cluster.idents=T,
        scale = F,
        group.by='seurat_clusters') +
  theme( axis.title.x = element_blank(),
         axis.title.y = element_text(size=20),
         axis.text.x = element_text(angle = 45, hjust = 1, size=12),
         axis.text.y = element_text(size=16),
         panel.spacing = unit(2, "lines"))+
  scale_y_discrete(limits=rev) +
  ylab('Cluster')
#===============================================================================



### Standard workflow, ging_sm
#===============================================================================
ging_sm <- NormalizeData(ging_sm)
ging_sm <- FindVariableFeatures(ging_sm)
ging_sm <- ScaleData(ging_sm, vars.to.regress = 'percent.mt')
ging_sm <- RunPCA(ging_sm, npcs=39)
ElbowPlot(ging_sm, ndims = 39)
ging_sm <- FindNeighbors(ging_sm,
                         reduction='pca',
                         dims=1:10)
ging_sm <- FindClusters(ging_sm,
                        resolution=0.5,
                        algorithm = 4)
ging_sm <- RunUMAP(ging_sm,
                   reduction='pca',
                   dims=1:10)
#===============================================================================

### Annotating clusters, ging_sm
#===============================================================================
umap_ging_sm <- DimPlot2(
  ging_sm,
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
umap_ging_sm

cd_ging_sm <- ClusterDistrBar(
  ging_sm$condition,
  ging_sm$subset_ann,
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
cd_ging_sm
#===============================================================================



### Standard workflow, ging_myeloid
#===============================================================================
ging_myeloid <- NormalizeData(ging_myeloid)
ging_myeloid <- FindVariableFeatures(ging_myeloid)
ging_myeloid <- ScaleData(ging_myeloid, vars.to.regress = 'percent.mt')
ging_myeloid <- RunPCA(ging_myeloid)
ElbowPlot(ging_myeloid, ndims = 50)
ging_myeloid <- FindNeighbors(ging_myeloid,
                              reduction='pca',
                              dims=1:10)
ging_myeloid <- FindClusters(ging_myeloid,
                             resolution=0.5,
                             algorithm = 4)
ging_myeloid <- RunUMAP(ging_myeloid,
                        reduction='pca',
                        dims=1:10)
#===============================================================================

### ORA, ging_myeloid
#===============================================================================
ging_myeloid <- JoinLayers(ging_myeloid)

ging_myeloid_de <- FindMarkers(ging_myeloid,
                               group.by = 'condition',
                               ident.1 = 'Vit D-',
                               ident.2 = 'Vit D+')

uni_genes_ging <- rownames(m_sob.list_mms_ging)
de_genes_ging_myeloid <- rownames(
  subset(ging_myeloid_de, p_val_adj <= 0.05)
)

ging_myeloid_ora <- enrichGO(de_genes_ging_myeloid,
                             OrgDb = org.Mm.eg.db,
                             universe = uni_genes_ging,
                             keyType = 'SYMBOL',
                             readable = T,
                             ont = 'ALL',
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.1)

ging_myeloid_ora_simp <- simplify(ging_myeloid_ora,
                                  cutoff=0.4,
                                  by='p.adjust',
                                  select_fun = min,
                                  measure='Wang')

clusterProfiler::dotplot(ging_myeloid_ora_simp,
                         showCategory = 10)
#===============================================================================


