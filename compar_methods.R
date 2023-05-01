## MS592 Data Analysis
# library(dplyr)
# library(Matrix)
# library(ComplexHeatmap)
# library(circlize)
# library(readxl)
# 
# setwd('~/BU/Sem2/MA592/Project/')
# 
# # read in seurat data
# srt_dat <- read_excel('KOvsControl_markers.xlsx') %>% as.data.frame() #read.csv('KOvsControl_markers.csv')
# 
# # read in FR-perturb data
# B <- read.csv('B_r5.csv', header = FALSE)
# B_p <- read.csv('pvals_r5.csv', header = FALSE)
# 
# p_names <- readLines('CD14_p_rownames.txt')
# g_names <- readLines('CD14_c_rownames.txt')
# 
# rownames(B) <- p_names
# colnames(B) <- g_names
# rownames(B_p) <- p_names
# colnames(B_p) <- g_names
# 
# # read in CI-IV data
# ci_dat <- read.csv('results_rep500.csv')

# plot 

# set target gene for comparison

target <- 'MAP3K7'
genes <- srt_dat$gene[grep(paste0(target, '_vs_Control'), srt_dat$compar)]

ci_vec <- ci_dat[which(ci_dat$XZ == target),]
rownames(ci_vec) <- ci_vec$Y

compar_eff <- data.frame(seurat_lfc = srt_dat$avg_log2FC[grep(paste0(target, '_vs_Control'), srt_dat$compar)], 
                         FR_perturb = B[target, genes] %>% t(),
                         CI_IV = ci_vec[genes, 'LATE_median'])
compar_eff <- as.matrix(compar_eff)*-1
colnames(compar_eff) <- c('Seurat', 'FactorizeRecover', 'CI_IV')
rownames(compar_eff) <- genes

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))
h <- Heatmap(compar_eff, col = col_fun, column_title = target, name = 'effect')
h

compar_p <- data.frame(seurat_lfc = srt_dat$p_val_adj[grep(paste0(target, '_vs_Control'), srt_dat$compar)], 
                         FR_perturb = B_p[target,genes] %>% t(),
                       CI_IV = ci_vec[genes, 'LATE_pval'])
compar_p <- as.matrix(compar_p)
rownames(compar_p) <- genes
colnames(compar_p) <- c('Seurat', 'FactorizeRecover', 'CI_IV')
# compar_p <- compar_p
col_fun = colorRamp2(c(0, 1), c("seagreen", "white"))
col_fun(seq(-3, 3))

h
Heatmap(compar_p, row_order = row_order(h), column_order = column_order(h), col = col_fun, column_title = 'p-values', name = ' ')

## check datasets
# check <- readLines('CD14_norm_c_rownames.txt')
# all(srt_dat$X %in% check)


## check specific
ci_dat <- read.csv('results_target_rep500.csv')

target <- 'TRAF6'
genes <- ci_dat$Y[grep(target, ci_dat$XZ)]

ci_vec <- ci_dat[which(ci_dat$XZ == target),]
rownames(ci_vec) <- ci_vec$Y

srt_vec <- srt_dat[grep(paste0(target, '_vs_Control'), srt_dat$compar),]
rownames(srt_vec) <- srt_vec$gene

compar_eff <- data.frame(seurat_lfc = srt_vec[genes, 'avg_log2FC'],
                         FR_perturb = B[target, genes] %>% t(),
                         CI_IV = -1*ci_vec[genes, 'LATE_median'])
compar_eff <- as.matrix(compar_eff)*-1
colnames(compar_eff) <- c('Seurat', 'FactorizeRecover', 'CI_IV')
rownames(compar_eff) <- genes

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))
h <- Heatmap(compar_eff, col = col_fun, column_title = target, name = 'effect')
h

compar_p <- data.frame(seurat_lfc = srt_vec[genes, 'p_val_adj'],
                       FR_perturb = B_p[target,genes] %>% t(),
                       CI_IV = ci_vec[genes, 'LATE_pval'])
compar_p <- as.matrix(compar_p)
rownames(compar_p) <- genes
colnames(compar_p) <- c('Seurat', 'FactorizeRecover', 'CI_IV')
# compar_p <- compar_p
col_fun = colorRamp2(c(0, 1), c("seagreen", "white"))
col_fun(seq(-3, 3))

h
Heatmap(compar_p, row_order = row_order(h), column_order = column_order(h), col = col_fun, column_title = 'p-values', name = ' ')
