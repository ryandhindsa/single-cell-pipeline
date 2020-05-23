################################################################################
# Seurat-based pipeline for integrated single-cell dat from mult. genotypes    #
# Author: Ryan Dhindsa                                                         #
################################################################################

# Imports ----------------------------------------------------------------------
# library(sctransform)
library(Seurat)
library(yaml)
library(cowplot)
library(data.table)
library(ggplot2)
theme_set(theme_cowplot())


# Globals ----------------------------------------------------------------------
config <- yaml.load_file("~/Dropbox/scRNA/projects/hnrnpu_hippocampus/config.yml")


# Functions --------------------------------------------------------------------
ReadData <- function(sample.matrix, min.cells, min.features, project.name) {
    # Reads in 10X Chromium count matrices and returns Seurat object
    #
    # Args:
    #   sample.matrix: a tab-delimited sample matrix containing paths to data
    #       must also include sample.name, genotype, gender, and species cols
    #
    #   min.cells: include features detected in at least this many cells
    #   min.features: include cells where at least this many features are 
    #                 detected
    #   project.name: name of project
    #
    # Returns:
    #   obj.list: a list of Seurat objects
    samples <- fread(sample.matrix)
    
    obj.list <- list()
    
    for (i in 1:nrow(samples)) {
        message("Reading in 10X data matrix: ", samples[i]$path)
        temp.dat <- Read10X(data.dir = samples[i]$path)
        obj <- CreateSeuratObject(counts = temp.dat, 
                                  min.cells = min.cells, 
                                  project = project.name)
        
        # add meta data to Seurat object
        obj$sample.name <- samples[i]$sample.name
        obj$genotype <- samples[i]$genotype
        obj$gender <- samples[i]$gender
        obj$species <- samples[i]$species
        obj <- CalcPctMito(obj)
        
        obj.list[[i]] <- obj
    }
    
    names(obj.list) <- samples$sample.name
    
    return(obj.list)
}


CalcPctMito <- function(obj) {
    # Calculates percentage of mitochondrial genes in object
    # Args:
    #   obj: a Seurat object with raw counts
    # Returns:
    #   obj: returns same object with mitochondrial %age stored in meta.data
    species <- unique(obj$species)
    
    if(species == "mouse") {
        mito.genes <- grep(pattern = "^mt-", 
                           x = rownames(x = obj), 
                           value = TRUE)
        
    } else if(species=="human") {
        mito.genes <- grep(pattern = "^MT-", 
                           x = rownames(x = obj), 
                           value = TRUE)
    } else {
        stop(paste(species), "is not a recognized species. 
             Please choose mouse or human.")
    }
    
    percent.mito <- Matrix::colSums(x = GetAssayData(obj, 
                                                     slot = 'counts')
                                    [mito.genes, ]) / 
        Matrix::colSums(x = GetAssayData(obj, slot = 'counts'))
    
    obj[['percent.mito']] <- percent.mito
    return(obj)
    }


VlnPlots <- function(obj.list, feature, max.cutoff=Inf, min.cutoff=-Inf) {
    # Takes in lost of Seurat objects and produces violin plot for given feature
    #
    # Args:
    #   obj.list: list of Seurat objects
    #   feature: feature to plot
    #
    # Returns:
    #   p: violin plot
    
    df <- data.frame(sample=factor(), feature=double())
    
    for (i in 1:length(obj.list)) {
        curr.df <- data.frame(obj.list[[i]]@meta.data$sample.name, 
                              obj.list[[i]]@meta.data[, feature])
        
        df <- rbind(df, curr.df)
    }
    
    colnames(df) <- c("sample", "feature")
    
    
    p <- ggplot(df, aes(x=sample, y=feature)) +
        geom_violin(scale = "width",
                    adjust = 1,
                    trim = TRUE,
                    fill = "#00BFC4") +
        ylab(feature) +
        geom_jitter(height = 0, size = .1, alpha=0.2) + 
        geom_hline(yintercept = max.cutoff, linetype=2) +
        geom_hline(yintercept = min.cutoff, linetype=2) +
        ggtitle(feature)
    
    return(p)
}


VisualizeQC <- function(obj.list, out.dir, max.pct.mito, 
                        min.n.genes, max.n.genes, max.umi) {
    # Creates and save violin plots for following QC features:
    #   n.genes, n.umi, and pct.mito
    #
    # Args:
    #   obj.list: list of Seurat objects
    #   out.dir: output directory where plots will be saved
    #   max.pct.mito: mitochondrial DNA %age cutoff
    #   min.n.genes: min number of genes 
    #   max.n.genes: max number of genes
    #   max.umi: max nUMI
    #
    # Returns:
    #   None
    qc.dir <- file.path(out.dir, "qc_plots")
    dir.create(qc.dir)
    
    pct.mito <- VlnPlots(all.obj, "percent.mito", 
                         max.cutoff = config$qc$max_pct_mito)
    
    fname <- file.path(qc.dir, "pct.mito.pdf")
    save_plot(fname, pct.mito, base_width = 6)
    
    n.genes <- VlnPlots(all.obj, "nFeature_RNA", 
                        min.cutoff = config$qc$min_n_genes, 
                        max.cutoff = config$qc$max_n_genes)
    
    fname <- file.path(qc.dir, "n.genes.pdf")
    save_plot(fname, n.genes, base_width = 6)
    
    n.umi <- VlnPlots(all.obj, "nCount_RNA", 
                      max.cutoff = config$qc$max_umi)
    fname <- file.path(qc.dir, "n.umi.pdf")
    save_plot(fname, n.umi, base_width = 6)
}


FilterData <- function(obj.list, min.n.genes, max.n.genes, 
                       max.pct.mito, max.umi) {
    # Filters cells in each seurat object contained in obj.list for nGenes, 
    #   nUMI, and pct mito
    # Args:
    #   obj.list: list of objects
    #   min.n.genes: min number of genes
    #   max.n.genes: max number of genes
    #   max.pct.mito: max percent mitochondrial dna
    # Returns:
    #   obj.list: filtered objects
    for (i in 1:length(obj.list)) {
        obj.list[[i]] <- subset(all.obj[[i]],
                                subset=nFeature_RNA > config$qc$min_n_genes &
                                    nFeature_RNA < config$qc$max_n_genes &
                                    percent.mito < config$qc$max_pct_mito &
                                    nCount_RNA < config$qc$max_umi)
    }
    
    return(obj.list)
}


NormalizeObj <- function(obj.list) {
    # Normalizes raw data in each object and finds varaible genes
    #  
    # Args:
    #   obj.list: a list of objects
    #
    # Returns:
    #   obj.list: list of objects with normalized data
    for (i in 1:length(x = obj.list)) {
        obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
        obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                              selection.method = "vst",
                                              nfeatures = 3000,
                                              verbose = FALSE)
    }
    return(obj.list)
}


IntegrateObj <- function(obj.list, n.dims) {
    # Identifies anchors and integrates normalized objects
    # Args:
    #   obj.list: list of normalized objects
    # Returns:
    #   integrated.obj: an integrated object
    anchors <- FindIntegrationAnchors(object.list = all.obj, dims = 1:n.dims)
    
    integrated.obj <- IntegrateData(anchorset = anchors, dims = 1:n.dims)
    DefaultAssay(object = integrated.obj) <- "integrated"
    
    # Run the standard workflow for visualization and clustering
    integrated.obj <- ScaleData(object = integrated.obj, 
                                vars.to.regress = c("nCount_RNA", 
                                                    "percent.mito"))
    integrated.obj <- RunPCA(object = integrated.obj, 
                             npcs = n.dims, verbose = FALSE)
    return(integrated.obj)
}


VisualizePCA <- function(integrated.obj, out.dir) {
    # Saves pdfs for pc loadings, heatmap, and elbow plot
    # Args: 
    #   integrated.obj: integrated seurat object
    #   out.dir: path to output directory
    # Returns:
    #   None
    pca.dir <- file.path(out.dir, "pca")
    dir.create(pca.dir)
    
    # pc loadings
    pc.loadings <- VizDimLoadings(integrated.obj, dims = 1:2)
    fname <- file.path(pca.dir, "pc_loadings.pdf")
    save_plot(fname, pc.loadings, base_width = 10, base_height = 7)
    
    # heatmap
    fname <- file.path(pca.dir, "pc_heatmap.pdf")
    pdf(fname, width = 8, height = 6)
    DimHeatmap(integrated.obj, dims = 1:6, cells = 500, 
                             balanced = TRUE)
    dev.off()
    
    # elbow plot
    pc.elbows <- ElbowPlot(object = integrated.obj, ndims = 30)
    fname <- file.path(pca.dir, "elbow_plot.pdf")
    save_plot(fname, pc.elbows, base_width = 6, base_height = 5)
}


ClusterIntegratedObj <- function(integrated.obj, n.dims, 
                                 clustering.res, umap.seed, save.obj=T) {
    # Clusters integrated data using SLM algorithm and performs UMAP for viz
    #   and saves objects as an rds object
    # Args:
    #   integrated.obj: integrated object
    #   n.dims: number of dimensions for clustering
    #   clustering.res: resolution for SLM
    #   umap.seed: seed for UMAP
    # Returns:
    #   integrated.obj: integrated object w/ cluster annotations
    integrated.obj <- FindNeighbors(object = integrated.obj, dims = 1:n.dims)
    integrated.obj <- FindClusters(object = integrated.obj, 
                                   resolution = clustering.res, algorithm = 3)
    
    integrated.obj <- RunUMAP(object = integrated.obj, reduction = "pca",
                   dims = 1:n.dims, seed.use = 20)

    if(save.obj) {
        saveRDS(integrated.obj, file = config$output$integrated_obj)
    }

    return(integrated.obj)
}


VisualizeClusters <- function(integrated.obj, out.dir) {
    # Saves pdfs for pc loadings, heatmap, and elbow plot
    # Args: 
    #   integrated.obj: integrated seurat object
    #   out.dir: path to output directory
    # Returns:
    #   None
    cluster.dir <- file.path(out.dir, "clustering")
    dir.create(cluster.dir)
    
    # pc loadings
    p1 <- DimPlot(object = obj, reduction = "umap", 
                  group.by = "genotype")
    
    fname <- file.path(cluster.dir, "genotype_umap.pdf")
    save_plot(fname, p1, base_width = 5.5, base_height = 4)
    
    # elbow plot
    p2 <- DimPlot(object = obj, reduction = "umap", 
                  group.by = "ident", 
                  label = TRUE, repel = T, pt.size = 0.01) + NoLegend()
    fname <- file.path(cluster.dir, "cluster_umap.pdf")
    save_plot(fname, p2, base_width = 5, base_height = 4)
}

# Main -------------------------------------------------------------------------
all.obj <- ReadData(config$sample_matrix, 
                    config$qc$min_cells, config$qc$min_features, 
                    config$project_name)

# Create violin plots
dir.create(config$output$dir)

VisualizeQC(obj.list = all.obj, out.dir = config$output$dir, 
            min.n.genes = config$qc$min_n_genes,
            max.n.genes = config$qc$max_n_genes,
            max.pct.mito = config$qc$max_pct_mito,
            max.umi = config$qc$max_umi)

# Filter
all.obj <- FilterData(all.obj,
                      min.n.genes = config$qc$min_n_genes,
                      max.n.genes = config$qc$max_n_genes,
                      max.pct.mito = config$qc$max_pct_mito,
                      max.umi = config$qc$max_umi)

# Normalize
all.obj <- NormalizeObj(all.obj)

# Integrate
obj <- IntegrateObj(all.obj, n.dims = config$integration$dims_to_calc)

VisualizePCA(obj, config$output$dir)

# Cluster
obj <- ClusterIntegratedObj(obj, 
                            n.dims = 11, 
                            clustering.res = 0.7, 
                            umap.seed = config$clustering$umap_seed, save.obj = F)

obj <- RunUMAP(object = obj, reduction = "pca", 
               dims = 1:11, seed.use = 1)

DimPlot(obj, label=T)

VisualizeClusters(integrated.obj = obj, out.dir = config$output$dir)


p1 <- DimPlot(object = obj, reduction = "umap", 
              group.by = "genotype")
p2 <- DimPlot(object = obj, reduction = "umap", 
              group.by = "ident", 
              label = TRUE, repel = T) + NoLegend()
plot_grid(p1, p2)
