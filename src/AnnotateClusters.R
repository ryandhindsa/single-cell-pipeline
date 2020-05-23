################################################################################
# Seurat-based pipeline for integrated single-cell dat from mult. genotypes    #
# Author: Ryan Dhindsa                                                         #
################################################################################

# Imports ----------------------------------------------------------------------
library(Seurat)
library(yaml)
library(cowplot)
library(data.table)
library(wesanderson)


# Globals ----------------------------------------------------------------------
config <- yaml.load_file("~/Dropbox/scRNA/projects/dnm1_p5/config.yml")
obj <- readRDS(config$output$integrated_obj)

# Functions --------------------------------------------------------------------
RenameClusterGroups <- function(obj, cluster.names.file, obj.path) {
    # Renames clusters
    # Args:
    #   obj: Seurat object
    #   cluster.names.file: path to txt file containing old and new cluster ids
    # Returns:
    #   obj: seurat object with new cluster identities
    
    # stash current identities
    obj[["orig.clusters"]] <- Idents(object = obj)
    
    cluster.names <- fread(cluster.names.file, 
                           header = T)
    
    current.cluster.ids <- cluster.names$old
    new.cluster.ids <- cluster.names$new

    names(new.cluster.ids) <- current.cluster.ids
    
    obj <- RenameIdents(obj, new.cluster.ids)
    
    # stash new identities in meta data
    obj[["celltype"]] <- obj@active.ident
    obj[["celltype.genotype"]] <- paste0(obj@active.ident, "_", obj$genotype)
    
    saveRDS(obj, obj.path)
    
    return(obj)
}


PlotClusters <- function(obj, out.dir) {
    # Re-plot UMAP w/ annotated clusters
    #
    # Args:
    #   obj: Seurat obj
    #   out.dir: output directory where plot will be saved
    #  
    # Returns:
    #   plot
    # out.dir <- config$output$dir
    cluster.dir <- file.path(out.dir, "clustering")
    suppressWarnings(dir.create(cluster.dir))
    
    p2 <- DimPlot(object = obj, reduction = "umap", 
                  group.by = "ident", 
                  label = T, repel = T) + NoLegend()
    fname <- file.path(cluster.dir, "annotated_cluster_umap.pdf")
    save_plot(fname, p2, base_width = 5, base_height = 4) 
    return(p2)
}

# Main -------------------------------------------------------------------------
cluster.names <- config$clustering$annotations
obj <- RenameClusterGroups(obj.integrated, cluster.names, config$output$integrated_obj)
PlotClusters(obj, config$output$dir)
DimPlot(obj, reduction = "umap", label = T, repel = T)

