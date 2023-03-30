library(Seurat)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

scanpy_to_seurat <- function(opref) {
    
    tbl = paste0(opref, '_x.tsv')
    metadata = paste0(opref, '_obs.tsv')

    # load normalized counts matrix
    X = read.table(tbl, header=TRUE) %>% t

    # load metadata table
    meta = read.table(metadata, header=TRUE)
    rownames(meta) <- meta$bc

    # create seurat object
    seurat_obj <- CreateSeuratObject(
        counts = X,
        data = X,
        meta = meta
    )
    Idents(seurat_obj) <- seurat_obj$sample


    # set factor level for leiden clusters:
    seurat_obj$leiden <- factor(seurat_obj$leiden, levels=unique(seurat_obj$leiden)[order(unique(seurat_obj$leiden))])

    umap <- as.matrix(select(seurat_obj@meta.data, c(umap_x, umap_y)))
    rownames(umap) <- colnames(seurat_obj)
    colnames(umap) = c('umap_1', 'umap_2')

    seurat_obj@reductions$umap <-  CreateDimReducObject(
        embeddings=umap,
        key="umap_",
        assay="RNA"
    )

    fname = paste0(opref, '.rds')
    saveRDS(seurat_obj, file=fname)
    return(seurat_obj)

}