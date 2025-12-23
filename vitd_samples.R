#### Introduction
#===============================================================================
# Title: mms_vitd_analysis.R
# Purpose: Import files, filter, add metadata
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

### Directories
#===============================================================================
# Directories created on machine for analysis. Count matrices were stored in 
# 'data' folder and this path was stored in 'data_path' variable.
figure_path <- here('figures')
anal_path <- here('analysis')
data_path <- here('data')
result_path <- here('results')
scripts_path <- here('scripts')

#===============================================================================

### Samples
#===============================================================================
buccalnt <- file.path(
  data_path,
  '20230822-B147S-MMC-GloverBuccalNT-1116HEF_20231218183536_TCM.tsv'
  )

buccaltx <- file.path(
  data_path,
  '20230822-B147S-MMC-GloverBuccalTx-1116HEF_20231218184133_TCM.tsv'
  )

gingivant <- file.path(
  data_path,
  '20230822-B147S-MMC-GloverGingivaNT-1116HEF_20231218182431_TCM.tsv'
  )

gingivatx <- file.path(
  data_path,
  '20230822-B147S-MMC-GloverGingivaTx-1116HEF_20231218182642_TCM.tsv'
  )

sample_list_tcm <- c(
  buccalnt,
  buccaltx,
  gingivant,
  gingivatx
  )
#===============================================================================

### Produce matrices
#===============================================================================
mat_list <- list()
sl <- sample_list_rcm

for (i in seq_along(sl)) {
  f <- sl[i]
  print(sl[i])
  m <- read.delim(f, 
                  header=TRUE, 
                  row.names=1, 
                  check.names=FALSE
  )
  m_name <- paste0('mat_', tools::file_path_sans_ext(basename(f)))
  mat_list[[m_name]]<- m
}
#===============================================================================

### Create Seurat objects
#===============================================================================
sob_list <- list()
sob_names <- list()

for (i in seq_along(mat_list)) {
  cur_mat <- mat_list[[i]]
  p_name <- gsub('mat_', '', names(mat_list[i]))
  p_name <- gsub('_RCM', '', p_name)
  print(p_name)
  sob <- CreateSeuratObject(counts = cur_mat,
                            project = p_name)
  sob_name <- paste0('sob_', p_name)
  sob_list[[sob_name]] <- sob
  assign(sob_name, sob)
}

sob_buccal_tx$orig.ident <- 'sob_buccal_tx'
sob_buccal_nt$orig.ident <- 'sob_buccal_nt'
sob_gingiva_tx$orig.ident <- 'sob_gingiva_tx'
sob_gingiva_nt$orig.ident <- 'sob_gingiva_nt'

#===============================================================================

### Annotation of metadata
#===============================================================================
# Vit D+ means sufficient Vit-D in diet
# Vit D- means insufficient Vit-D in diet
sob_buccal_nt$condition <- 'Vit D-'
sob_buccal_tx$condition <- 'Vit D+'
sob_gingiva_nt$condition <- 'Vit D-'
sob_gingiva_tx$condition <- 'Vit D+'
#===============================================================================

### Combine Seurat objects 
#===============================================================================
# making combined 1) buccal, 2) gingiva

# Create buccal obj
m_sob_mms_buccal <- merge(
  sob_buccal_tx,
  y = sob_buccal_nt
  )

# Create gingival obj
m_sob_mms_ging <- merge(
  sob_gingiva_tx,
  y = sob_gingiva_nt
  )

#===============================================================================

### Filtering
#===============================================================================
m_sob_mms_buccal$percent.mt<- PercentageFeatureSet(
  m_sob_mms_buccal,
  pattern='^mt-'
)

m_sob_mms_buccal_f <- subset(
  m_sob_mms_buccal,
  subset = nFeature_RNA>=50 &
    nFeature_RNA<=2500 &
    percent.mt<=15
)

m_sob.list_mms_buccal <- m_sob_mms_buccal_f




m_sob_mms_ging$percent.mt<- PercentageFeatureSet(
  m_sob_mms_ging,
  pattern='^mt-'
)

m_sob_mms_ging_f <- subset(
  m_sob_mms_ging,
  subset = nFeature_RNA>=50 &
    nFeature_RNA<=2500 &
    percent.mt<=20
)

m_sob.list_mms_ging <- m_sob_mms_ging_f
#===============================================================================

### Annotation of merged objects
#===============================================================================
m_sob.list_mms_buccal@meta.data$condition[m_sob.list_mms_buccal@meta.data$orig.ident == 'sob_buccal_nt'] <- 'Vit D-'
m_sob.list_mms_buccal@meta.data$condition[m_sob.list_mms_buccal@meta.data$orig.ident == 'sob_buccal_tx'] <- 'Vit D+'

m_sob.list_mms_ging@meta.data$condition[m_sob.list_mms_ging@meta.data$orig.ident == 'sob_gingiva_nt'] <- 'Vit D-'
m_sob.list_mms_ging@meta.data$condition[m_sob.list_mms_ging@meta.data$orig.ident == 'sob_gingiva_tx'] <- 'Vit D+'
#===============================================================================
