# Loading the libraries.
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dowser))
suppressPackageStartupMessages(library(scoper))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidyverse))

#**********************************|Importing the data from a repository|*********************************#
# Importing the data from a local repository.
Merged_bcr <- readChangeoDb('Path_to_directory')
Merged_bcr <- Merged_bcr %>% filter(productive)

# Importing GEX data as a seurat object for integrating with VDJ data
gex_b_sub_so_5_mito <- readRDS('Path_to_directory')
#*********************************************************************************************************#

#******************************|Filtering the dataset for downstream analysis|****************************#
# Removing cells with multiple heavy chain.
multi_heavy <- table(filter(Merged_bcr, locus == "IGH")$cell_id)
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]
Merged_bcr <- filter(Merged_bcr, !cell_id %in% multi_heavy_cells)

# Spliting cells by heavy and light chains
heavy_cells <- filter(Merged_bcr, locus == "IGH")$cell_id
light_cells <- filter(Merged_bcr, locus == "IGK" | locus == "IGL")$cell_id
no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]
Merged_bcr <- filter(Merged_bcr, !cell_id %in% no_heavy_cells)

# Inspect the data object
# The gene expression Seurat object.
print(gex_b_sub_so_5_mito)

# Cell type annotations
head(Idents(gex_b_sub_so_5_mito),1)

# The immune repertoire data.
head(Merged_bcr,1)
#*********************************************************************************************************#

#***************************************|Merging VDJ and GEX datasets|************************************#
# merging information from GEX dataset with vdj dataset.
#Changing cell ids for meta data table
md <- gex_b_sub_so_5_mito@meta.data
rownames(md) = gsub("TP1_", "", rownames(md))
rownames(md) = gsub("TP2_", "", rownames(md))
rownames(md) = gsub("TP3_", "", rownames(md))
rownames(md) = gsub("TP4_", "", rownames(md))
rownames(md) = gsub("TP5_", "", rownames(md))
head(md)
gex_b_sub_so_5_mito@meta.data = md

#Creating matching index again for removal of cells that does not matches.
match.index <- match(Merged_bcr$cell_id, rownames(gex_b_sub_so_5_mito@meta.data))

#Estimating the proportion of the cells do not matched.
mean(is.na(match.index))

#Transferring the gex information into the vdj dataset.
cells.sample = as.character(gex_b_sub_so_5_mito$orig.ident)
Merged_bcr$sample = unlist(lapply(match.index, function(x){ifelse(is.na(x), NA, cells.sample[x])}))
cells.annotations = as.character(gex_b_sub_so_5_mito$custom.cell.types)
Merged_bcr$annotations = unlist(lapply(match.index, function(x){ifelse(is.na(x), NA, cells.annotations[x])}))
cells.clusters = as.character(gex_b_sub_so_5_mito$seurat_clusters)
Merged_bcr$seurat_clusters = unlist(lapply(match.index, function(x){ifelse(is.na(x), NA, cells.clusters[x])}))

# Making a new index column for matching with a dataset.
Merged_bcr <- Merged_bcr %>% rowwise() %>%
  mutate(cell_id_unique = paste(sample, cell_id, sep = "_"))

# Removing the transformed seurat object and re uploading.
rm(gex_b_sub_so_5_mito)
gex_b_sub_so_5_mito <- readRDS('Path_to_directory')

# check out select columns.
head(select(Merged_bcr, cell_id, v_call, j_call, sample, annotations, seurat_clusters),5)

# Make cell IDs in BCR match those in Seurat Object.
Merged_bcr$cell_id_unique = paste0(Merged_bcr$sample, "_", Merged_bcr$cell_id)
Merged_bcr$cell_id_unique[1]
#*********************************************************************************************************#

#**********************************|Separating heavy and light chains data|*******************************#
# Separating the light and heavy chains from the dataset.
Merged_bcr_L <- subset(Merged_bcr, Merged_bcr$locus == "IGL")
Merged_bcr_K <- subset(Merged_bcr, Merged_bcr$locus == "IGK")
Merged_bcr_L <- rbind(Merged_bcr_L, Merged_bcr_K)

Merged_bcr_H <- subset(Merged_bcr, Merged_bcr$locus == "IGH")
#*********************************************************************************************************#

#**********************************|Step wise Calculation of SHM distance|********************************#
# Creating a v_length and j_length column for light chains in the dataset.
Merged_bcr_L$v_length <- (Merged_bcr_L$v_sequence_end - Merged_bcr_L$v_sequence_start) + 1
Merged_bcr_L$j_length <- (Merged_bcr_L$j_sequence_end - Merged_bcr_L$j_sequence_start) + 1
Merged_bcr_L$v_dist <- 1 - Merged_bcr_L$v_identity
Merged_bcr_L$j_dist <- 1 - Merged_bcr_L$j_identity
Merged_bcr_L$v_shm_dist <- Merged_bcr_L$v_length * Merged_bcr_L$v_dist 
Merged_bcr_L$j_shm_dist <- Merged_bcr_L$j_length * Merged_bcr_L$j_dist 
Merged_bcr_L$shm_values <- Merged_bcr_L$v_shm_dist + Merged_bcr_L$j_shm_dist

# Creating a v_length and j_length column for heavy chains in the dataset.
Merged_bcr_H$v_length <- (Merged_bcr_H$v_sequence_end - Merged_bcr_H$v_sequence_start) + 1
Merged_bcr_H$j_length <- (Merged_bcr_H$j_sequence_end - Merged_bcr_H$j_sequence_start) + 1
Merged_bcr_H$v_dist <- 1 - Merged_bcr_H$v_identity
Merged_bcr_H$j_dist <- 1 - Merged_bcr_H$j_identity
Merged_bcr_H$v_shm_dist <- Merged_bcr_H$v_length * Merged_bcr_H$v_dist 
Merged_bcr_H$j_shm_dist <- Merged_bcr_H$j_length * Merged_bcr_H$j_dist 
Merged_bcr_H$shm_values <- Merged_bcr_H$v_shm_dist + Merged_bcr_H$j_shm_dist
#*********************************************************************************************************#

#***********************************|Adding labels and removing NA values|********************************#
# First removing rows that has NA values in seurat column.
Merged_bcr_L <- Merged_bcr_L[!is.na(Merged_bcr_L$seurat_clusters),]
Merged_bcr_H <- Merged_bcr_H[!is.na(Merged_bcr_H$seurat_clusters),]

# Adding labels for the tree tips.
labels <- read.delim("Path_to_directory/tips_labels.txt", header = T)

Merged_bcr_H <- merge(Merged_bcr_H, labels, by = "cell_id", all.x = TRUE)
Merged_bcr_H  <- Merged_bcr_H [!is.na(Merged_bcr_H$tip_labels),]
#*********************************************************************************************************#

#****************************************|Identifying clonal clusters|************************************#
# Picking a threshold using shazam.
# Finding threshold using heavy chains
Merged_bcr_H <- filter(Merged_bcr, locus=="IGH")
dist_ham <- distToNearest(Merged_bcr_H)
head(dist_ham) %>%
  select(cell_id_unique, dist_nearest)
output <- findThreshold(dist_ham$dist_nearest)
threshold <- output@threshold

# Visualizing the distance-to-nearest distribution and threshold
plotDensityThreshold(output)

# Performing clustering using scoper
# Assigning the clonal clusters
results <- hierarchicalClones(dist_ham, threshold=threshold)
results_db <- as.data.frame(results)
head(results_db) %>%
  select(cell_id_unique, clone_id)

# Visualizing clone size distribution.
# Getting clone sizes using dplyr functions
clone_sizes <- countClones(results_db)

# Plotting cells per clone
ggplot(clone_sizes, aes(x=seq_count))+
  geom_bar() + theme_bw() +
  xlab("Sequences per clone")
#*********************************************************************************************************#

#*******************************|Reconstructing clonal germlines using dowser|****************************#
#Reading IMGT reference file from a local folder.
references <- readIMGT(dir = "Path_to_directory")

#Reconstructing germlines.
results_db = createGermlines(results_db, references)

#Checking the output.
results_db$germline_alignment_d_mask[1]
#*********************************************************************************************************#

#**************************************|Building and visualizing tree|************************************#
# Formatting clones with dowser. Making clone objects with aligned, processed sequences and collapsing 
# identical sequences unless differ by trait. Adding up duplicate_count column for collapsed sequences
# store day, isotype, gex_annotation. Discarding clones with < 5 distinct sequences. 
clones = formatClones(results_db,
                      traits = c("sample", "c_call", "annotations", "cell_id_unique", "seurat_clusters"),
                      minseq=5)
clones

#Collecting cell ids for top three clonotypes.
ids_clons1 <- clones$data[[1]]
ids_clons1 <- as.data.frame(ids_clons1@data)
ids_clons1 <- ids_clons1[!is.na(ids_clons1$sample),]
ids_clons1$clonotype <- "clonotype_1"

ids_clons2 <- clones$data[[2]]
ids_clons2 <- as.data.frame(ids_clons2@data) 
ids_clons2 <- ids_clons2[!is.na(ids_clons2$sample),]
ids_clons2$clonotype <- "clonotype_2"

ids_clons3 <- clones$data[[3]]
ids_clons3 <- as.data.frame(ids_clons3@data) 
ids_clons3 <- ids_clons3[!is.na(ids_clons3$sample),]
ids_clons3$clonotype <- "Clonotype_3"

top_3_clons <- rbind(ids_clons1, ids_clons2, ids_clons3)

# Using maximum likelihood method.
trees = getTrees(clones)
head(trees)

# Plotting all trees
plots = plotTrees(trees, tips = "sample", tipsize=2)
plots[[1]]

p1 <- plotTrees(trees, scale = 0.1)[[1]] +
  geom_tiplab(aes(label = tip_labels)) +
  geom_tippoint(aes(fill=factor(c_call)), size = 2.5, color = "black", shape = 21) + 
  scale_fill_manual(values = c('#E24832', '#FDDB8B', '#578BBF'))

p2 <- plotTrees(trees, scale = 0.1)[[2]] +
  geom_tippoint(aes(fill=factor(c_call)), size = 2.5, color = "black", shape = 21) +
  scale_fill_manual(values = c("#D93629", "#F99355", "#FEE699", "#E6F5EC", "#95C7DF", "#497AB6"))

p3 <- plotTrees(trees, scale = 0.1)[[3]] +
  geom_tippoint(aes(fill=factor(c_call)), size = 2.5, color = "black", shape = 21) +
  scale_fill_manual(values = c("#D93629", "#F99355", "#FEE699", "#E6F5EC", "#95C7DF", "#497AB6"))
#*********************************************************************************************************#
















