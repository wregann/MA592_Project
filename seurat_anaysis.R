library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(remotes)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ComplexHeatmap)
library(Matrix)

setwd("~/BU/Cleary/CompressTest")

# install_github("mojaveazure/seurat-disk")

dat <- readRDS('../Data/DougPaper/Normal_seurat.rds')
dat@meta.data$Cell <- rownames(dat@meta.data)

setwd('~/BU/Sem2/MA592/Project/')

c_rownames <- readLines('CD14_c_rownames.txt')
c_colnames <- readLines('CD14_c_colnames.txt')

p_rownames <- readLines('CD14_p_rownames.txt')
p_colnames <- readLines('CD14_p_colnames.txt')

c <- dat@assays$RNA@data[c_rownames, c_colnames]
p <- dat@assays$perturbations[p_rownames,p_colnames]
meta <- dat@meta.data[c_colnames,] 

# cells to subset to

cells <- dat@meta.data$Cell
drop_cells <- c()
keep_cells <- c()
p_full <- dat@assays$perturbations

for (i in 1:ncol(p_full)){
  perts <- rownames(p_full)[which(p_full[,i] > 0)]
  
  if (all(perts %in% p_rownames)){
    keep_cells <- c(keep_cells, colnames(p_full)[i])
  }
}

# create new object

CD14.dat <- CreateSeuratObject(
  c,
  project = "CD14_Inf",
  assay = "RNA",
  meta.data = meta
)
CD14.dat[["perturbations"]] <- CreateAssayObject(counts = p)
# CD14.dat <- subset(x = dat, subset = Cell %in% c_colnames)
# CD14.dat@assays$perturbations <- CD14.dat@assays$perturbations[p_rownames,p_colnames]
saveRDS(CD14.dat, file = 'CD14_dat.rds')

### write files
setwd('norm/')
writeMM(CD14.dat@assays$RNA@data, file = 'CD14_norm_counts.mtx')
writeLines(CD14.dat@assays$RNA@data %>% rownames(), file('CD14_norm_c_rownames.txt'))
writeLines(CD14.dat@assays$RNA@data %>% colnames(), file('CD14_norm_c_colnames.txt'))

writeMM(CD14.dat@assays$perturbations@counts, file = 'CD14_norm_p.mtx')
writeLines(CD14.dat@assays$perturbations@counts %>% rownames(), file('CD14_norm_p_rownames.txt'))
writeLines(CD14.dat@assays$perturbations@counts %>% colnames(), file('CD14_norm_p_colnames.txt'))

##
setwd('~/BU/Sem2/MA592/Project/')
CD14.dat <- readRDS('CD14_dat.rds')
# write.csv(CD14.dat@meta.data, 'CD14_metadata.csv')
mat <- CD14.dat@assays$RNA@counts %>% as.matrix()
# mat <- mat[c_rownames,]
dat <- CD14.dat@assays$RNA@data
P <- CD14.dat@assays$perturbations

gene_cor <- cor(mat %>% t())
Heatmap(gene_cor)

# object2 <- CreateSeuratObject(mat) # Create a new Seurat object with just the genes of interest
# orig.ident <- CD14.dat@meta.data # Pull the identities from the original Seurat object as a data.frame
# object2 <- AddMetaData(object = object2, metadata = orig.ident) # Add the idents to the meta.data slot
# object2@assays$perturbations <- p

all.genes <- rownames(CD14.dat@assays$RNA@data)
CD14.dat <- ScaleData(CD14.dat, features = all.genes)
CD14.dat <- FindVariableFeatures(object = CD14.dat)
CD14.dat <- RunPCA(CD14.dat, features = VariableFeatures(object = CD14.dat))
CD14.dat <- FindNeighbors(CD14.dat, dims = 1:10)
CD14.dat <- FindClusters(CD14.dat, resolution = 0.5)
CD14.dat <- RunUMAP(CD14.dat, dims = 1:10)
DimPlot(CD14.dat, reduction = "umap")
FeaturePlot(CD14.dat, features = 'n.guides')

group <- CD14.dat@meta.data$guide.names
group <- group %>% 
  gsub(pattern = 'safe-targeting', replacement = 'Control', .) %>%
  gsub(pattern = 'non-targeting', replacement = 'Control', .) %>%
  gsub(pattern = 'Control--', replacement = '', .) %>%
  gsub(pattern = '--Control', replacement = '', .)
# p <- CD14.dat@assays$perturbations
# group <- rep('NA', ncol(p))
# for (i in 1:ncol(p)){
#   perts <- rownames(p)[which(p[,i] > 0)]
#   perts <- perts[-which(perts %in% c("non-targeting", "safe-targeting"))]
#   g <- paste(perts, collapse = '_')
#   group[i] <- g
# }
# group[which(group == '')] <- 'Control'
CD14.dat@meta.data$Group <- group

group_list <- group %>% unique() %>% sort()
col_pal <- rep('white', 92)
col_pal[which(group_list == 'Control')] <- 'slategray2'
col_pal[which(group_list == 'TRAF6')] <- 'steelblue2'
col_pal[which(group_list == 'IKBKB')] <- 'springgreen2'
col_pal[which(group_list == 'NFKB1')] <- 'tan'
col_pal[which(group_list == 'MAPK1')] <- 'purple2'
col_pal[which(group_list == 'MAP3K7')] <- 'pink1'

library(RColorBrewer)
my_cols = brewer.pal(8,"Dark2")
tplot <- DimPlot(main.integrated, group.by='group',cols=alpha(my_cols,0.66),pt.size=1)
ko_dim <- DimPlot(CD14.dat, reduction = "umap", group.by = 'Group', cols = alpha(col_pal, 0.6), pt.size = 4)
        # cells.highlight = which(group %in% c('Control', 'TRAF6', 'IKBKB', 'NFKB1', 'MAPK1', 'MAP3K7')),
        # sizes.highlight = 4)
ko_dim
umap_df <- CD14.dat[["umap"]]@cell.embeddings %>% as.data.frame()
umap_df$id <- rownames(umap_df)
meta_df <- CD14.dat@meta.data
meta_df$id <- rownames(meta_df)
merge_df <- merge(meta_df, umap_df, by = 'id')

ggplot(merge_df[-which(merge_df$Group == 'Control'),], aes(x=UMAP_1, y=UMAP_2, col=Group)) +
  geom_point(size=1, alpha = 0.5)

markers <- FindAllMarkers(CD14.dat)
write.csv(markers, file = 'merkers_res0.5.csv')

markers <- FindAllMarkers(CD14.dat)
write.csv(markers, file = 'KO_markers.csv')

FeaturePlot(CD14.dat, features = c('CD14'))

df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df) <- c('p_val', 'avg_log2FC', 'pct.1', 'pct.2',    'p_val_adj', 'compar', 'gene')
Idents(CD14.dat) <- 'Group'
compar <- unique(group)
compar <- compar[-which(compar == 'Control')]
for (i in compar){
  if (length(which(CD14.dat@meta.data$Group == i)) < 5){next}
  print(i)
  g.markers <- FindMarkers(CD14.dat, ident.1 = i, ident.2 = "Control", logfc.threshold = 0.05)
  g.markers <- g.markers[which(g.markers$p_val_adj < 1),]
  if (nrow(g.markers) == 0){next}
  g.markers$compar <- paste0(i, '_vs_Control')
  g.markers$gene <- rownames(g.markers)
  df <- rbind(df, g.markers)
}
all(rownames(df) %in% rownames(CD14.dat@assays$RNA@data))
all(df$gene %in% rownames(CD14.dat@assays$RNA@data))
write.csv(df, file = 'KOvsControl_markers.csv')

g.markers <- FindMarkers(CD14.dat, ident.1 = "NFKBIA", ident.2 = "Control")
# SetIdent(CD14.dat, 'Group')
