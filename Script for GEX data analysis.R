# Loading libraries that are required for the analysis.
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(djvdj))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(celldex))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tmaptools))
suppressPackageStartupMessages(library(gghalves))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(dittoSeq))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(monocle3))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(scCustomize))

# Setting an additional function for making table containing statistics of cells.
averagenCountnFeature <- function(cells) {
  output <- tibble(
    'sample' = character(),
    'mean(nCount)' = numeric(),
    'median(nCount)' = numeric(),
    'mean(nFeature)' = numeric(),
    'median(nFeature)' = numeric(),
    'mean(percent_mito)' = numeric()
  )
  for ( i in levels(cells$sample) ) {
    tmp <- tibble(
      'sample' = i,
      'mean(nCount)' = cells %>% dplyr::filter(sample == i) %>% pull(nCount) %>% mean(),
      'median(nCount)' = cells %>% dplyr::filter(sample == i) %>% pull(nCount) %>% median(),
      'mean(nFeature)' = cells %>% dplyr::filter(sample == i) %>% pull(nFeature) %>% mean(),
      'median(nFeature)' = cells %>% dplyr::filter(sample == i) %>% pull(nFeature) %>% median(),
      'mean(percent_mito)' = cells %>% dplyr::filter(sample == i) %>% pull(percent_mito) %>% mean()
    )
    output <- bind_rows(output, tmp)
  }
  return(output)
}

# Making customized color palletes.
pal_npg_4 = pal_npg("nrc")(4)
pal_npg_5 = pal_npg("nrc")(5)
pal_npg_6 = pal_npg("nrc")(6)
pal_npg_7 = pal_npg("nrc")(7)
pal_npg_10 = pal_npg("nrc")(10)

#**************************|Importing local files and creating seurat objects|***************************#
# Creating a Seurat object that contain only RNA assay. 
#Specifying directories containing gene expression matrices acquired from cellranger processing.
TP1_rawdata <- Read10X(data.dir = "path_to_directory")
TP2_rawdata <- Read10X(data.dir = "path_to_directory")
TP3_rawdata <- Read10X(data.dir = "path_to_directory")
TP4_rawdata <- Read10X(data.dir = "path_to_directory")
TP5_rawdata <- Read10X(data.dir = "path_to_directory")

#Creating Seurat object from the raw data for gene expression.
TP1 <- CreateSeuratObject(counts = TP1_rawdata$`Gene Expression`, project = "TP1")
TP2 <- CreateSeuratObject(counts = TP2_rawdata$`Gene Expression`, project = "TP2")
TP3 <- CreateSeuratObject(counts = TP3_rawdata$`Gene Expression`, project = "TP3")
TP4 <- CreateSeuratObject(counts = TP4_rawdata$`Gene Expression`, project = "TP4")
TP5 <- CreateSeuratObject(counts = TP5_rawdata$`Gene Expression`, project = "TP5")
#********************************************************************************************************#

#***********************|Filtering and processing cells for downstream analysis|*************************#
# Merging all time points together.
gex_so <- merge(TP1, y = c(TP2, TP3, TP4, TP5), add.cell.ids = c("TP1", "TP2", "TP3", "TP4", "TP5"), project = "Matze")
gex_so

# Making a cell features counting table.
#Counting mitochondrial gene expression.
gex_so[["percent.mt"]] <- PercentageFeatureSet(gex_so, pattern = "^MT-")
cells <- tibble(
  cell = colnames(gex_so),
  sample = as.factor(gex_so$orig.ident),
  nCount = gex_so$nCount_RNA,
  nFeature = gex_so$nFeature_RNA,
  percent_mito = gex_so$percent.mt
)
raw.cells.table <- averagenCountnFeature(cells) %>% knitr::kable()
raw.cells.table

# Plotting a violin plot for QC metrics before filtering.
p1 <- ggplot(cells, aes(x = sample, y = nCount, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_point(position = position_jitter(width = 0.15), aes(colour = sample)) +
  theme_bw() +
  scale_fill_manual(values = pal_npg_5) +
  scale_color_manual(values = pal_npg_5) +
  scale_x_discrete(limits = levels(cells$sample)) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = 'Transcripts expression') +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )

p2 <- ggplot(cells, aes(x = sample, y = nFeature, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_point(position = position_jitter(width = 0.15), aes(colour = sample)) +
  theme_bw() +
  scale_fill_manual(values = pal_npg_5) +
  scale_color_manual(values = pal_npg_5) +
  scale_x_discrete(limits = levels(cells$sample)) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = 'Unique genes expression') +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )

p3 <- ggplot(cells, aes(x = sample, y = percent_mito, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_point(position = position_jitter(width = 0.15), aes(colour = sample)) +
  theme_bw() +
  scale_fill_manual(values = pal_npg_5) +
  scale_color_manual(values = pal_npg_5) +
  scale_x_discrete(limits = levels(cells$sample)) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = 'Mitochondrial genes expression') +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )

# Making a features counting table for cells
#applying a cutoff on mitochondrial gene expression (<5%).
gex_so_5_mito <- subset(gex_so, subset = percent.mt < 5)
cells <- tibble(
  cell = colnames(gex_so_5_mito),
  sample = as.factor(gex_so_5_mito$orig.ident),
  nCount = gex_so_5_mito$nCount_RNA,
  nFeature = gex_so_5_mito$nFeature_RNA,
  percent_mito = gex_so_5_mito$percent.mt
)
filt.cells.table <- averagenCountnFeature(cells) %>% knitr::kable()
filt.cells.table

# Plotting a violin plot for QC metrics after filtering.
p1 <- ggplot(cells, aes(x = sample, y = nCount, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_point(position = position_jitter(width = 0.15), aes(colour = sample)) +
  theme_bw() +
  scale_fill_manual(values = pal_npg_5) +
  scale_color_manual(values = pal_npg_5) +
  scale_x_discrete(limits = levels(cells$sample)) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = 'Transcripts expression') +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )

p2 <- ggplot(cells, aes(x = sample, y = nFeature, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_point(position = position_jitter(width = 0.15), aes(colour = sample)) +
  theme_bw() +
  scale_fill_manual(values = pal_npg_5) +
  scale_color_manual(values = pal_npg_5) +
  scale_x_discrete(limits = levels(cells$sample)) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = 'Unique genes expression') +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )

p3 <- ggplot(cells, aes(x = sample, y = percent_mito, fill = sample)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_point(position = position_jitter(width = 0.15), aes(colour = sample)) +
  theme_bw() +
  scale_fill_manual(values = pal_npg_5) +
  scale_color_manual(values = pal_npg_5) +
  scale_x_discrete(limits = levels(cells$sample)) +
  scale_y_continuous(labels = scales::comma) +
  labs(y = 'Mitochondrial genes expression') +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  )

# Normalization, finding the most variable genes and scaling the data for downstrema analysis. 
#Normalization.
gex_so_5_mito <- NormalizeData(gex_so_5_mito, normalization.method = "LogNormalize", scale.factor = 10000)

#Finding variable genes.
gex_so_5_mito <- FindVariableFeatures(gex_so_5_mito, selection.method = "vst")
#********************************************************************************************************#

#*************************|Applying dimensional reduction methods and plotting|**************************#
#Scaling the count data.
gex_so_5_mito <- ScaleData(gex_so_5_mito, features = rownames(gex_so_5_mito))

# Running PCA dimensionality reduction method and calculation of PCs.
#Running PCA.
gex_so_5_mito <- RunPCA(gex_so_5_mito, features = rownames(gex_so_5_mito), npcs = 50)

#Plotting elbowplot for visualizing StDv.
ElbowPlot(gex_so_5_mito)

#In addition, also calculate by taking a K-mean of 10. (So, taking 15 PCs is an ideal number.)
intrinsicDimension::maxLikGlobalDimEst(gex_so_5_mito@reductions$pca@cell.embeddings, k = 10)

# Using the computed number of PCs for finding nearest clusters and then further running UMAP and tSNEs.
#Finding neigbors and clustring the cells.
gex_so_5_mito <- FindNeighbors(gex_so_5_mito, dims = 1:15)
gex_so_5_mito <- FindClusters(gex_so_5_mito, resolution = 0.7)

#Running UMAP.
gex_so_5_mito <- RunUMAP(gex_so_5_mito, dims = 1:15)

#Running tSNE.
gex_so_5_mito <- RunTSNE(gex_so_5_mito, dims = 1:15)

# Visualizing clusters on UMAP labelled with cluster number and timepoints.
UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = gex_so_5_mito@meta.data$seurat_clusters
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))
p1 <- bind_cols(gex_so_5_mito@meta.data, as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.5) +
  geom_label(
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 3.0,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.25,
    label.size = 0,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = pal_npg_10,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(gex_so_5_mito@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  ) +
  labs(title = "UMAP for cell clusters")

UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = gex_so_5_mito@meta.data$orig.ident
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))
p2 <- bind_cols(gex_so_5_mito@meta.data, as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = orig.ident)) +
  geom_point(size = 0.4) +
  geom_label(
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 3.5,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.25,
    label.size = 0,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = pal_npg_5,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(gex_so_5_mito@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  ) +
  labs(title = "UMAP for timepoints")
#********************************************************************************************************#

#*********************|Assignment of cell type using Human Primary Cell Atlas Data|**********************#
# Annotating cell types based on a single cell expression.
#Retrieving Human Primary Cell Atlas dataset as a reference.
singler_ref <- HumanPrimaryCellAtlasData()

#Extracting signature of expression for main labels.
singler_results_hpca_main_5 <- SingleR::SingleR(
  test = GetAssayData(gex_so_5_mito, assay = 'RNA', slot = 'data'),
  ref = singler_ref,
  labels = singler_ref@colData@listData$label.main
)

p <- plotScoreHeatmap(
  singler_results_hpca_main_5,
  show.labels = TRUE,
  annotation_col = data.frame(
    donor = gex_so_5_mito@meta.data$orig.ident,
    row.names = rownames(singler_results_hpca_main_5)
  )
)

# Assignment of subcell types basaed on distinguished cell type expression.
singler_results_hpca_fine_5 <- SingleR::SingleR(
  test = GetAssayData(gex_so_5_mito, assay = 'RNA', slot = 'data'),
  ref = singler_ref,
  labels = singler_ref@colData@listData$label.fine
)

p <- plotScoreHeatmap(
  singler_results_hpca_fine_5,
  show.labels = F,
  annotation_col = data.frame(
    donor = gex_so_5_mito@meta.data$orig.ident,
    row.names = rownames(singler_results_hpca_fine_5)
  )
)

# Storing labels and scores in the seurat object.
gex_so_5_mito@meta.data$cell_type_singler_hpca_main_5 <- singler_results_hpca_main_5@listData$labels
gex_so_5_mito@meta.data$cell_type_singler_hpca_fine_5 <- singler_results_hpca_fine_5@listData$labels
singler_scores_main_5 <- singler_results_hpca_main_5@listData$scores %>%
  as_tibble() %>%
  dplyr::mutate(assigned_score = NA)
singler_scores_fine_5 <- singler_results_hpca_fine_5@listData$scores %>%
  as_tibble() %>%
  dplyr::mutate(assigned_score = NA)
for (i in seq_len(nrow(singler_scores_main_5))) {
  singler_scores_main_5$assigned_score[i] <- singler_scores_main_5[[singler_results_hpca_main_5@listData$labels[i]]][i]
}
for (i in seq_len(nrow(singler_scores_fine_5))) {
  singler_scores_fine_5$assigned_score[i] <- singler_scores_fine_5[[singler_results_hpca_fine_5@listData$labels[i]]][i]
}
gex_so_5_mito@meta.data$cell_type_singler_hpca_main_score_5 <- singler_scores_main_5$assigned_score
gex_so_5_mito@meta.data$cell_type_singler_hpca_fine_score_5 <- singler_scores_fine_5$assigned_score

# Plotting UMAP to visualize the cell labels at main level.
UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = gex_so_5_mito@meta.data$cell_type_singler_hpca_main_5
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

p1 <- bind_cols(gex_so_5_mito@meta.data, as.data.frame(gex_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = cell_type_singler_hpca_main_5)) +
  geom_point(size = 0.5) +
  geom_label(
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 3.5,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.3,
    label.size = 0,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = pal_npg_4,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(gex_so_5_mito@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  ) +
  labs(title = "UMAP for cell types")

# Making a table for cell types annotated by SingleR method.
table_samples_by_cell_type_5 <- gex_so_5_mito@meta.data %>%
  dplyr::group_by(orig.ident, cell_type_singler_hpca_main_5) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(cell_type_singler_hpca_main_5, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("orig.ident", "total_cell_count", dplyr::everything()))

table_samples_by_cell_type_5 %>% knitr::kable()
table_samples_by_cell_type_5

# Plotting the composition of the clusters as stacked barplots.
temp_labels <- gex_so_5_mito@meta.data %>%
  group_by(orig.ident) %>%
  tally()

p1 <- table_samples_by_cell_type_5 %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(cluster = factor(orig.ident, levels = levels(gex_so_5_mito@meta.data$orig.ident))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity', colour = "black") +
  geom_text(
    data = temp_labels,
    aes(
      x = orig.ident,
      y = Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ), color = 'black', size = 1.5
  ) +
  scale_fill_manual(name = 'Cell type', values = pal_npg_4) +
  scale_y_continuous(
    name = 'Number of cells',
    labels = scales::comma,
    expand = c(0.01,0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  ) + labs(title = "Composition of cell types")

# Removing cells other than B cells that are annotated by SingleR at main and fine level.
gex_b_so_5_mito <- subset(gex_so_5_mito, subset = cell_type_singler_hpca_main_5 == "B_cell")
Idents(gex_b_so_5_mito) <- "cell_type_singler_hpca_fine_5"
gex_b_sub_so_5_mito <- subset(gex_b_so_5_mito, idents = c("NK_cell:CD56hiCD62L+", "T_cell:CD4+_central_memory"), invert = T)

# Plotting a UMAP for visualizing cells annotated at fine (subtype) level.
UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = gex_b_sub_so_5_mito@meta.data$cell_type_singler_hpca_fine_5
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

p1 <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = cell_type_singler_hpca_fine_5)) +
  geom_point(size = 0.5) +
  geom_label(
    data = UMAP_centers_cluster,
    mapping = aes(x, y, label = cluster),
    size = 3.5,
    fill = 'white',
    color = 'black',
    fontface = 'bold',
    alpha = 0.3,
    label.size = 0,
    show.legend = FALSE
  ) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = pal_npg_7,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right') +
  coord_fixed() +
  annotate(
    geom = 'text', x = Inf, y = -Inf,
    label = paste0('n = ', format(nrow(gex_b_sub_so_5_mito@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  ) +
  labs(title = "UMAP for B cells subtypes")

# Making a table for automatically assigned B cell subtypes via SingleR.
table_samples_by_cell_type_5 <- gex_b_sub_so_5_mito@meta.data %>%
  dplyr::group_by(orig.ident, cell_type_singler_hpca_fine_5) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(cell_type_singler_hpca_fine_5, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("orig.ident", "total_cell_count", dplyr::everything()))

table_samples_by_cell_type_5 %>% knitr::kable()
table_samples_by_cell_type_5

# Plotting the composition as stacked barplots.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(orig.ident) %>%
  tally()

p1 <- table_samples_by_cell_type_5 %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(cluster = factor(orig.ident, levels = levels(gex_b_sub_so_5_mito@meta.data$orig.ident))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity', colour = "black") +
  geom_text(
    data = temp_labels,
    aes(
      x = orig.ident,
      y = Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ), color = 'black', size = 1.5
  ) +
  scale_fill_manual(name = 'Cell type', values = pal_npg_7) +
  scale_y_continuous(
    name = 'Number of cells',
    labels = scales::comma,
    expand = c(0.01,0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  ) + labs(title = "Composition of B cell subtypes")
#********************************************************************************************************#

#**************************************|Second round of clustering|**************************************#
# Running a second round of clustering after removing non-B cells from the samples.
#Finding variable genes.
gex_b_sub_so_5_mito <- FindVariableFeatures(gex_b_sub_so_5_mito, selection.method = "vst")

#Scaling the count data.
gex_b_sub_so_5_mito <- ScaleData(gex_b_sub_so_5_mito, features = rownames(gex_b_sub_so_5_mito))

#Scaling the count data.
gex_b_sub_so_5_mito <- ScaleData(gex_b_sub_so_5_mito, features = rownames(gex_b_sub_so_5_mito))

#Running PCA.
gex_b_sub_so_5_mito <- RunPCA(gex_b_sub_so_5_mito, features = rownames(gex_b_sub_so_5_mito), npcs = 50)

#Plotting elbowplot for visualizing StDv.
ElbowPlot(gex_b_sub_so_5_mito)

#Finding neighbors and clusters.
gex_b_sub_so_5_mito <- FindNeighbors(gex_b_sub_so_5_mito, dims = 1:10)
gex_b_sub_so_5_mito <- FindClusters(gex_b_sub_so_5_mito, resolution = 0.20)

#Running UMAP.
gex_b_sub_so_5_mito <- RunUMAP(gex_b_sub_so_5_mito, dims = 1:10)
#********************************************************************************************************#

#*****************|Identifying differentially expressed genes and plotting volcano plots|****************#
#Making violin plots.
#Finding markers in specific cluster.
markers_plasmavsBcells <- FindMarkers(gex_b_sub_so_5_mito, ident.1 = c(2, 5), ident.2 = c(0, 1, 3, 4), logfc.threshold = 0.25, test.use = "DESeq2")
markers_MBCsBVac_vs_otherBcells <- FindMarkers(gex_b_sub_so_5_mito, ident.1 = 4, ident.2 = c(0, 1, 3), logfc.threshold = 0.25, test.use = "DESeq2")
markers_MBCsAVac_vs_otherBcells <- FindMarkers(gex_b_sub_so_5_mito, ident.1 = 0, ident.2 = c(1, 3, 4), logfc.threshold = 0.25, test.use = "DESeq2")
markers_MBCsCXAVac_vs_otherBcells <- FindMarkers(gex_b_sub_so_5_mito, ident.1 = 3, ident.2 = c(0, 1, 4), logfc.threshold = 0.25, test.use = "DESeq2")
markers_Atypical_vs_otherBcells <- FindMarkers(gex_b_sub_so_5_mito, ident.1 = 1, ident.2 = c(0, 3, 4), logfc.threshold = 0.25, test.use = "DESeq2")

#Building Volcano plots.
p1 <- EnhancedVolcano(markers_plasmavsBcells, lab = rownames(markers_plasmavsBcells), x = 'avg_log2FC', y = 'p_val_adj', 
                      FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                      col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                      legendLabels=c('Not sig.','log2FC','adj.P',
                                     'adj.P & log2FC'), border = 'full', borderWidth = 0.5, 
                      labCol = '#000000', selectLab = c("CD38", "MS4A1", "MKI67", "XBP1", "CXCR3"), boxedLabels = TRUE, colAlpha = 4/5, labFace = 'bold',
                      legendLabSize = 10, labSize = 4, xlim = c(-5,5), ylim = c(0,300), 
                      title = "Plasma vs B cells", subtitle = NULL, colConnectors = 'black', drawConnectors = TRUE,widthConnectors = 1.0) + theme(aspect.ratio = 1)

p2 <- EnhancedVolcano(markers_MBCsBVac_vs_otherBcells, lab = rownames(markers_MBCsBVac_vs_otherBcells), x = 'avg_log2FC', y = 'p_val_adj', 
                      FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                      col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                      legendLabels=c('Not sig.','log2FC','adj.P',
                                     'adj.P & log2FC'), border = 'full', borderWidth = 0.5, 
                      labCol = '#000000', selectLab = c("IRF1", "GRASP", "NFKBID", "ACTG1", "ACTB"), boxedLabels = TRUE, colAlpha = 4/5, labFace = 'bold',
                      legendLabSize = 10, labSize = 4, xlim = c(-4,4), ylim = c(0,100), 
                      title = "MBCs before vac. vs other B cells", subtitle = NULL, colConnectors = 'black', drawConnectors = TRUE,widthConnectors = 1.0) + theme(aspect.ratio = 1)

p3 <- EnhancedVolcano(markers_MBCsAVac_vs_otherBcells, lab = rownames(markers_MBCsAVac_vs_otherBcells), x = 'avg_log2FC', y = 'p_val_adj', 
                      FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                      col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                      legendLabels=c('Not sig.','log2FC','adj.P',
                                     'adj.P & log2FC'), border = 'full', borderWidth = 0.5, 
                      labCol = '#000000', selectLab = c("NFKBIA", "LTB", "ACTG1"), boxedLabels = TRUE, colAlpha = 4/5, labFace = 'bold',
                      legendLabSize = 10, labSize = 4, xlim = c(-4,4), ylim = c(0,75), 
                      title = "MBCs after vac. vs other B cells", subtitle = NULL, colConnectors = 'black', drawConnectors = TRUE,widthConnectors = 1.0) + theme(aspect.ratio = 1)

p4 <- EnhancedVolcano(markers_MBCsCXAVac_vs_otherBcells, lab = rownames(markers_MBCsCXAVac_vs_otherBcells), x = 'avg_log2FC', y = 'p_val_adj', 
                      FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                      col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                      legendLabels=c('Not sig.','log2FC','adj.P',
                                     'adj.P & log2FC'), border = 'full', borderWidth = 0.5, 
                      labCol = '#000000', selectLab = c("ACTG1", "S100A4", "SELL", "CXCR3", "TXNIP"), boxedLabels = TRUE, colAlpha = 4/5, labFace = 'bold',
                      legendLabSize = 10, labSize = 4, xlim = c(-4,4), ylim = c(0,200), 
                      title = "MBCs after vac. CXCR3+ vs other B cells", subtitle = NULL, colConnectors = 'black', drawConnectors = TRUE,widthConnectors = 1.0) + theme(aspect.ratio = 1)

p5 <- EnhancedVolcano(markers_Atypical_vs_otherBcells, lab = rownames(markers_Atypical_vs_otherBcells), x = 'avg_log2FC', y = 'p_val_adj', 
                      FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                      col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                      legendLabels=c('Not sig.','log2FC','adj.P',
                                     'adj.P & log2FC'), border = 'full', borderWidth = 0.5, 
                      labCol = '#000000', selectLab = c("FCRL5", "FCRL3", "FGR", "ITGAX", "NFKBIA"), boxedLabels = TRUE, colAlpha = 4/5, labFace = 'bold',
                      legendLabSize = 10, labSize = 4, xlim = c(-4,4), ylim = c(0,180), 
                      title = "Atypical B cells vs other B cells", subtitle = NULL, colConnectors = 'black', drawConnectors = TRUE,widthConnectors = 1.0) + theme(aspect.ratio = 1)
#********************************************************************************************************#

#*************|Annotating cell clusters on the basis of marker genes expression and plotting|************#
# Observing marker gene expression in the clusters
your_markers <- c("CD19", "CD79A", "CD79B", "CD38", "CXCR3", "FCRL5", "FCRL3")
p <- dittoDotPlot(gex_b_sub_so_5_mito, your_markers, group.by = "seurat_clusters", main = "Marker gene expression across the clusters", min.color = "#1E90FF", max.color = "#DC143C",
                  min = -1, max = 1) + coord_flip()

# Building a data frame for customized cluster numbers and assigned cell type names.
clusters <- c(0, 1, 2, 3, 4, 5) 
custom.cell.types <- c("MBCs after vacc.", "Atypical BCs", "Plasma cells", "CXCR+ MBCs after vacc.", "MBCs before vacc.", "Plasma cells")
clust_no <- c('0', '1', '2', '3', '4', '5')
clusters.rename <- data.frame(clusters, custom.cell.types, clust_no)
gex_b_sub_so_5_mito@meta.data$custom.cell.types <- clusters.rename[gex_b_sub_so_5_mito@meta.data$seurat_clusters,]$custom.cell.types

# Creating a UMAP for visualizing custom cell type annotations.
UMAP_centers_cluster <- tibble(
  UMAP_1 = as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)$UMAP_1,
  UMAP_2 = as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)$UMAP_2,
  cluster = gex_b_sub_so_5_mito@meta.data$custom.cell.types
) %>%
  group_by(cluster) %>%
  summarize(x = median(UMAP_1), y = median(UMAP_2))

p1 <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = custom.cell.types)) +
  geom_point(size = 1.5) +
  theme_bw() +
  scale_color_manual(
    name = 'Cluster', values = pal_npg_5,
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  theme(legend.position = 'right', aspect.ratio = 1) +
  coord_fixed() +
  annotate(
    geom = 'text', x = -Inf, y = Inf,
    label = paste0('n = ', format(nrow(gex_b_sub_so_5_mito@meta.data), big.mark = ',', trim = TRUE)),
    vjust = -1.5, hjust = 1.25, color = 'black', size = 2.5
  ) +
  labs(title = "UMAP of cell types")

# Creating a table containing number of B subtypes.
table_samples_by_cell_type_5 <- gex_b_sub_so_5_mito@meta.data %>%
  dplyr::group_by(orig.ident, custom.cell.types) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::spread(custom.cell.types, count, fill = 0) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c("orig.ident", "total_cell_count", dplyr::everything()))
table_samples_by_cell_type_5 %>% knitr::kable()
table_samples_by_cell_type_5

# Visualizing a stacked barplot containing composition of B cell subtypes at each timepoint.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(orig.ident) %>%
  tally()

p1 <- table_samples_by_cell_type_5 %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  mutate(cluster = factor(orig.ident, levels = levels(gex_b_sub_so_5_mito@meta.data$orig.ident))) %>%
  ggplot(aes(orig.ident, value)) +
  geom_bar(aes(fill = variable), position = 'stack', stat = 'identity', colour = "black") +
  geom_text(
    data = temp_labels,
    aes(
      x = orig.ident,
      y = Inf,
      label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)),
      vjust = -1
    ), color = 'black', size = 1
  ) +
  scale_fill_manual(name = 'Cell type', values = pal_npg_5) +
  scale_y_continuous(
    name = 'Number of cells',
    labels = scales::comma,
    expand = c(0.01,0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  ) + labs(title = "Composition of custom annoations of B cell subtypes.")

#Highlighting cells in UMAP by sample
#Get cell names
Idents(gex_b_sub_so_5_mito) <- "orig.ident"
TP1 <- WhichCells(gex_b_sub_so_5_mito, idents = 'TP1')
TP2 <- WhichCells(gex_b_sub_so_5_mito, idents = 'TP2')
TP3 <- WhichCells(gex_b_sub_so_5_mito, idents = 'TP3')
TP4 <- WhichCells(gex_b_sub_so_5_mito, idents = 'TP4')
TP5 <- WhichCells(gex_b_sub_so_5_mito, idents = 'TP5')

#Make into list
cells_TP1 <- list(TP1 = TP1)
cells_TP2 <- list(TP2 = TP2)
cells_TP3 <- list(TP3 = TP3)
cells_TP4 <- list(TP4 = TP4)
cells_TP5 <- list(TP5 = TP5)

#Plotting as UMAPs
p1 <- Cell_Highlight_Plot(gex_b_sub_so_5_mito, cells_highlight = cells_TP1, highlight_color = "red") + labs(title = "Timepoint 1")
p2 <- Cell_Highlight_Plot(gex_b_sub_so_5_mito, cells_highlight = cells_TP2, highlight_color = "red") + labs(title = "Timepoint 2")
p3 <- Cell_Highlight_Plot(gex_b_sub_so_5_mito, cells_highlight = cells_TP3, highlight_color = "red") + labs(title = "Timepoint 3")
p4 <- Cell_Highlight_Plot(gex_b_sub_so_5_mito, cells_highlight = cells_TP4, highlight_color = "red") + labs(title = "Timepoint 4")
p5 <- Cell_Highlight_Plot(gex_b_sub_so_5_mito, cells_highlight = cells_TP5, highlight_color = "red") + labs(title = "Timepoint 5")

# Plotting UMAPs to visualize features 
#For all timepoints.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(orig.ident) %>%
  tally()
all_timepoints <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = orig.ident)) +
  geom_point(size = 0.4, show.legend = F) +
  theme_bw() +
  scale_color_manual(values = pal_npg_5) +
  labs(color = 'Timepoint') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  theme(aspect.ratio = 1,
        strip.text = element_text(face = 'bold', size = 10)
  )  + labs(title = "UMAPs for Timepoints")

#By timepoints as a split.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(orig.ident) %>%
  tally()
timepoints <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = orig.ident)) +
  geom_point(size = 0.4, show.legend = T) +
  theme_bw() +
  scale_color_manual(values = pal_npg_5) +
  labs(color = 'Timepoint') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  facet_wrap(~orig.ident, ncol = 5) +
  theme( aspect.ratio = 1,
         legend.position = 'right',
         strip.text = element_text(face = 'bold', size = 10)
  )

#For all Clusters.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()
all_clusters <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.4, show.legend = F) +
  theme_bw() +
  scale_color_manual(values = pal_npg_6) +
  labs(color = 'Clusters') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  theme(aspect.ratio = 1,
        strip.text = element_text(face = 'bold', size = 10)
  ) + labs(title = "UMAPs for Timepoints with clusters")

#Clusters at timepoints as a split.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()
clusters <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = seurat_clusters)) +
  geom_point(size = 0.4, show.legend = T) +
  theme_bw() +
  scale_color_manual(values = pal_npg_6) +
  labs(color = 'Clusters') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  facet_wrap(~orig.ident, ncol = 5) +
  theme(aspect.ratio = 1,
        legend.position = 'right',
        strip.text = element_text(face = 'bold', size = 10)
  )

# For all Celltypes at timepoints as a split.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(custom.cell.types) %>%
  tally()
all_cells <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = custom.cell.types)) +
  geom_point(size = 0.4, show.legend = F) +
  theme_bw() +
  scale_color_manual(values = pal_npg_5) +
  labs(color = 'Cell_types') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  theme(aspect.ratio = 1,
        strip.text = element_text(face = 'bold', size = 10)
  ) + labs(title = "UMAPs for Timepoints with ccell types")


#Celltypes at timepoints as a split.
temp_labels <- gex_b_sub_so_5_mito@meta.data %>%
  group_by(custom.cell.types) %>%
  tally()
cells <- bind_cols(gex_b_sub_so_5_mito@meta.data, as.data.frame(gex_b_sub_so_5_mito@reductions$umap@cell.embeddings)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color = custom.cell.types)) +
  geom_point(size = 0.4, show.legend = T) +
  theme_bw() +
  scale_color_manual(values = pal_npg_5) +
  labs(color = 'Cell_types') +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  facet_wrap(~orig.ident, ncol = 5) +
  theme(aspect.ratio = 1,
        legend.position = 'right',
        strip.text = element_text(face = 'bold', size = 10)
  )
#*******************************************************************************************************#

#************************************|Performing trajectory analysis|***********************************#
# Conversion of seurat object into cell dataset object.
cds <- as.cell_data_set(gex_b_sub_so_5_mito)
cds

# Getting cell metadata.
colData(cds)

# Observing gene names.
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_names <- rownames(fData(cds))

# Clustering cells (using clustering info from seurat's UMAP).
# Assigning partitions. 
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters$UMAP$partitions <- recreate.partition

# Assigning the cluster info.
list_cluster <- gex_b_sub_so_5_mito@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assigning UMAP coordiante or cell embeddings.
cds@int_colData@listData$reducedDims$UMAP <- gex_b_sub_so_5_mito@reductions$umap@cell.embeddings

# Applying trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

# Setting cluster 4 as an inception point.
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 4]))

# Visualizing as a plot
p1 <- plot_cells(cds,
                 color_cells_by = 'pseudotime',
                 label_groups_by_cluster = T,
                 label_branch_points = T,
                 label_roots = T,
                 label_leaves = T,
                 cell_size = 1.5,
                 graph_label_size = 3,
                 trajectory_graph_segment_size = 1.5) + theme(aspect.ratio = 1)
#*******************************************************************************************************#
