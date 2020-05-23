# Imports ----------------------------------------------------------------------
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(gprofiler2)
library(biomaRt)
theme_set(theme_cowplot())

# Objects ----------------------------------------------------------------------
setwd("~/Dropbox/IGM/hnrnpu_manuscript/")
hippo <- readRDS("results/hippocampus_subset.rds")
ctx <- readRDS("results/cortex_subset.rds")

# Global variables -------------------------------------------------------------
hippo.pal <- c("CA1" = "#EC579AFF", "CA2/CA3" = "#EE0011FF", 
               "Dentate Gyrus" = "#15983DFF", 
               "Subiculum" = "#0C5BB0FF", "Entorhinal Cortex" = "#FA6B09FF", 
               "IP" = "#149BEDFF", "Radial Glia" = "#A1C720FF", "SST" = "#E69F00", 
               "VIP" = "#16A08CFF", 
               "Prolif." = "black", "OPC"= "lightgrey", "Astrocyte"="red", "Cajal" = "blue")

ctx.pal <- c("Inhib. progenitor" = "#EC579AFF", "L2-4" = "#EE0011FF", 
             "L5/6" = "#15983DFF", "LGE" = "#0C5BB0FF", "SPN" = "#FA6B09FF", 
             "IP" = "#149BEDFF", "Radial Glia" = "#A1C720FF", "SST" = "#E69F00", 
             "VIP" = "#16A08CFF")


pal <- c(
  "Subiculum" = "#4E79A7", "Astrocyte" = "#79706E", "OPC" = "#D7B5A6",
  "IP" = "#A0CBE8", "Radial Glia" = "#FFBE7D", "Inhib. progenitor" = "#BAB0AC",
  "L2-4" = "#59A14F", "L5/6" = "#F28E2B", "Cajal" = "#499894", 
  "SPN" = "#D4A6C8", "SST" = "#E15759", "Entorhinal Cortex" = "#FF9D9A", 
  "LGE" = "#4E79A7", 
  "Dentate Gyrus" = "#B07AA1", "Prolif." = "#86BCB6", 
  "VIP" = "#F1CE63", "CA1" = "#59A14F", "CA2/CA3" = "#F28E2B")

FONT.SIZE <- 8
theme_set(theme_cowplot(font_size=FONT.SIZE)) # set font size for plots


# Functions --------------------------------------------------------------------
VisualizeUMAP <- function(obj) {
  # Function to create labeled UMAP without axes or legend
  #
  # Args:
  #   obj: a Seurat object
  #
  # Returns:
  #   p: a scatter plot of the first 2 UMAP dimensions
  df <- FetchData(obj, c("UMAP_1", "UMAP_2", "ident", "genotype", "sample.name"))
  
  p <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, fill=ident, col=ident)) + 
    geom_point(shape=16, size=.1) +
    scale_fill_manual(values=pal) +
    scale_color_manual(values=pal) + 
    xlab("UMAP 1") + 
    ylab("UMAP 2") + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) + 
    theme(legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=1)))
  
  p
}


CreateHeatmap <- function(obj, markers.path) {
  # Plot heatmap demonstrating expression of canonical markers
  # 
  # Args:
  #   obj: seurat object
  #   markers.path: path to text file containing markers (with colname "marker")
  #
  # Returns:
  #   p: a heatmap
  markers <- fread(markers.path)

  p <- DoHeatmap(subset(obj, downsample = 200), 
            features = markers$marker, 
            raster = F,
            size = 4) +
    NoLegend()
  
  p
}


CreateDotPlot <- function(obj, markers.path) {
  # Plot heatmap demonstrating expression of canonical markers
  # 
  # Args:
  #   obj: seurat object
  #   markers.path: path to text file containing markers (with colname "marker")
  #
  # Returns:
  #   p: a heatmap
  markers <- fread(markers.path)
  
  p <- DotPlot(obj, features = markers$marker, col.min = 0.5, 
               assay = "integrated", dot.scale = 4) + 
    RotatedAxis() +
    # coord_flip() +
    theme(text = element_text(size=FONT.SIZE), 
          axis.text=element_text(size=FONT.SIZE)) + 
    theme(legend.key.size = unit(.3, "cm")) + 
    guides(size = guide_legend(title = '% of cells')) + 
    guides(color = guide_colorbar(title = 'Expression')) +
    xlab("") + 
    ylab("") 
  
  p
}


'%!in%' <- function(x,y){
  # Function for declaring "is not"
  !('%in%'(x,y))
}
  

PerformDGE <- function(obj, celltype, min.pct = 0.1, logfc.threshold=0, 
                       genes=NULL, max.cells=Inf, test.use="MAST") {
  # Performs MAST DGE
  # 
  # Args:
  #   obj: seurat object
  #   celltype: celltype of interest
  #   min.pct: only test genes expressed in at least this proportion of cells
  #   logfc.threshold: test genes above this logfc
  #   genes: df of genes you want to include in DGE testing
  #
  # Returns:
  #   dge.results: df of DGE results
  DefaultAssay(obj) <- "RNA"
  obj$gender_binary <- ifelse(obj$gender == "M", 1, 0)
  
  obj$ngeneson <- scale(obj$nFeature_RNA)
  
  if(test.use == "MAST") {
    dge.results <- FindMarkers(obj, ident.1 = "HET", group.by = "genotype", 
                               subset.ident = celltype, assay = "RNA", 
                               features = genes, test.use = "MAST", 
                               logfc.threshold = logfc.threshold, 
                               min.pct = min.pct, 
                               max.cells.per.ident = max.cells,
                               latent.vars = c("ngeneson", 
                                               "gender_binary"))
    
  } else if (test.use == "MAST_downsample") {  
    dge.results <- FindMarkers(obj, ident.1 = "HET", group.by = "genotype", 
                               subset.ident = celltype, assay = "RNA", 
                               features = genes, test.use = "MAST", 
                               logfc.threshold = logfc.threshold, 
                               min.pct = min.pct, 
                               max.cells.per.ident = max.cells,
                               latent.vars = c("ngeneson"))
  } else {
    dge.results <- FindMarkers(obj, ident.1 = "HET", group.by = "genotype", 
                               subset.ident = celltype, assay = "RNA", 
                               features = genes,
                               logfc.threshold = logfc.threshold, 
                               min.pct = min.pct, 
                               max.cells.per.ident = max.cells)
  }
  
  dge.results <- dge.results %>% rownames_to_column('gene')
  
  dge.results$avg_log2FC <- dge.results$avg_logFC * log2(exp(1)) 
  
  dge.results$FDR <- p.adjust(dge.results$p_val, method = "BH")
  
  return(dge.results)
}


MultiDGE <- function(obj, min.pct=0.1, logfc.threshold=0, out.file, genes=NULL,
                     max.cells=Inf, exclude.celltypes=NULL, test.use="MAST") {
  # Performs DGE for each cell type contained in arg celltypes and writes
  # results to txt file
  #
  # Args:
  #   obj: Seurat object
  #   min.pct: only include genes expressed in at least this proportion of 
  #       cells
  #   logfc.threshold: only test genes w/ this log fc threshold
  #   out.dir: output directory
  #   genes: df of genes to include for DGE
  #
  # Returns:
  #   None
  
  all.dge <- list()
  for(i in unique(Idents(obj))){
    if(i %in% exclude.celltypes) {
      next
    } 
    message(paste("Performing DGE for celltype:", i))

    temp.dge <- PerformDGE(obj, i, min.pct, logfc.threshold, genes, 
                           max.cells = max.cells, test.use = test.use)
    
    # Write to list of dataframes
    temp.dge$celltype <- i
    all.dge[[i]] <- temp.dge
  }
  
  all.dge.df <- rbindlist(all.dge)
  all.dge.df$FDR <- p.adjust(all.dge.df$p_val, method="BH")
  
  write.table(all.dge.df, out.file, sep="\t", quote=F, row.names = F)
  
  return(all.dge.df)
}


BulkDGE <- function(obj, genes, min.pct = 0.1, logfc.threshold = 0) {
  # Run a pseudo-bulk RNA-seq analysis using MAST
  #
  # Args:
  #   obj: Seurat object
  #   genes: genes to include in the DGE analysis
  #   min.pct: test genes expressed in at least this  pct of cells
  #   
  # Returns:
  #   dge.results: bulk DGE results
  DefaultAssay(obj) <- "RNA"
  obj$gender_binary <- ifelse(obj$gender == "M", 1, 0)
  
  obj$ngeneson <- scale(obj$nFeature_RNA)
  
  dge.results <- FindMarkers(obj, ident.1 = "HET", group.by = "genotype", 
                             assay = "RNA", 
                             features = genes, test.use = "MAST", 
                             logfc.threshold = logfc.threshold, 
                             min.pct = min.pct, 
                             latent.vars = c("ngeneson", 
                                             "gender_binary"))
  
  dge.results <- dge.results %>% rownames_to_column('gene')
  dge.results$avg_log2FC <- dge.results$avg_logFC * log2(exp(1)) 
  dge.results$FDR <- p.adjust(dge.results$p_val, method = "BH")
  
  return(dge.results)
}


CountCells <- function(obj, group.by) {
  # Count the number of cells per celltype in each sample
  # 
  # Args:
  #   obj: seurat object
  #   group.by: group by (e.g. sample name or genotype)
  #
  # Returns:
  #   n.cells: dataframe containing the number of cells per celltype in each
  #     sample
  n.cells <- data.frame(table(Idents(obj), obj@meta.data[, group.by]))
  colnames(n.cells) <- c("Celltype", "Sample", "Count")
  return(n.cells)
}


PlotCellCounts <- function(obj, group.by = "sample.name") {
  # Create barchart of number of cells per identity
  #
  # Args:
  #   obj: Seurat object
  #   group.by: what to group counts by (e.g. genotype vs. sample name)
  # 
  # Returns:
  #   p: a barchart
  cell.counts <- CountCells(obj, group.by = group.by)
  p <- ggplot(cell.counts, aes(x=reorder(Celltype, Count), 
                               y=Count, fill=Sample)) + 
    geom_bar(stat="identity", position = "dodge") +
    scale_fill_brewer(palette="Paired") + 
    RotatedAxis() +
    xlab("")
  p
}


VolcanoPlot <- function(dge.df, log2fc.cutoff=0.14, fdr.cutoff=.05,
                        label.log2fc, label.fdr, color.pal) {
  # Create a volcano plot of DEGs
  #
  # Args:
  #   dge.df: a dataframe containing genes and FDRs
  #   log2fc.cutoff: cutoff for what you consider differentially expressed
  #     defaults to 10% change in expr.
  #   fdr.cutoff: FDR cutoff for DEG
  #   label.log2fc: log2fc cutoff for adding labels
  #   label.fdr: fdr cutoff for labels
  #   color.pal: a color palette 
  #
  # Returns:
  #   p: a volcano plot (ggplot2 object)
  
  dge.df <- subset(dge.df, abs(avg_log2FC) > log2fc.cutoff &
                     FDR < fdr.cutoff)
  
  p <- ggplot(dge.df, aes(x=avg_log2FC, y=-log10(FDR), color = celltype)) + 
    geom_vline(xintercept = 0, col="gray") +
    geom_hline(yintercept=-log10(.05), col="gray", linetype = "longdash") +
    geom_jitter(size=1, alpha=0.5) +
    geom_label_repel(aes(label=gene), size=2.5,
                     xlim = c(log2fc.cutoff, NA),
                     data=subset(dge.df, avg_log2FC > label.log2fc & 
                                   dge.df$FDR < label.fdr), 
                     show.legend = F, alpha=1) +
    geom_label_repel(aes(label=gene), size=2.5,
                     xlim = c(NA, -log2fc.cutoff),
                     data=subset(dge.df, avg_log2FC < -label.log2fc & 
                                   dge.df$FDR < label.fdr), 
                     show.legend = F, alpha=1) +
    scale_color_manual(values=color.pal) + 
    guides(color=guide_legend(ncol=2, 
                              override.aes = list(alpha=1))) +
    theme(legend.justification = c(1, 1), 
          legend.position = c(1, 1),
          legend.title = element_blank()) +
    xlab("log2(fold change)")
  
  p
}


PlotNumDEG <- function(dge.df, log2fc.cutoff=0.14, fdr.cutoff=0.05, color.pal) {
  # Plot number of DEGs per cell type
  #
  # Args:
  #   dge.df: dataframe of differentially expressed genes
  #   log2f.cutoff: log fc cutoff
  #   fdr.cutoff: fdr cutoff
  #   color.pal: color palette for plot
  #
  # Returns:
  #   p: ggplot object

  dge.tallies <- dge.df %>% 
    group_by(celltype) %>% 
    summarise(n=sum(abs(avg_log2FC) > log2fc.cutoff & FDR < fdr.cutoff))
  
  p <- ggplot(dge.tallies, aes(x=reorder(celltype, n), y=n, fill=celltype)) + 
    geom_bar(stat="identity", col="black", size=0.25) +
    scale_fill_manual(values=color.pal) + 
    NoLegend() +
    xlab("") +
    ylab("# of DEGs")
   
  p
}


SplitVlnPlot <- function(obj, gene.name, exclude.celltypes=NULL) {
  # Plot violin plots split by genotype
  # 
  # Args:
  #   obj: Seurat object
  #   gene.name: gene name
  #
  # Returns:
  #   p: a Seurat object
  DefaultAssay(obj) <- "RNA"
  
  dat <- FetchData(obj, 
                   vars = c(gene.name, "ident", "genotype", "sample.name"))
  
  dat <- dat %>% filter(ident %!in% exclude.celltypes)
  
  colnames(dat) <- c("Expression", "ident", "genotype", "sample.name")
  
  dat$genotype <- factor(dat$genotype, levels=c("WT", "HET"))
  
  p <- ggplot(dat, aes(x=reorder(ident, Expression), y=Expression, 
                     fill=genotype)) + 
    geom_split_violin(adjust = 1, scale = "width", size = 0.1) + 
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, 
                 geom = "errorbar",
                 width = 0.5,
                 position =  position_dodge(width = .5), 
                 col = "black") + 
    RotatedAxis() +
    scale_fill_manual(values = c("#86BCB6", "#FFBE7D")) + 
    xlab("") +
    ylab("Expression level") + 
    ggtitle(gene.name) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.title=element_blank())
  
  p 
}


JitterPlot <- function(obj, gene.name, group.by) {
  # Plot expression of a gene as a violin plot with overlaid points
  # 
  # Args:
  #   obj: Seurat object
  #   gene.name: gene name
  #   group.by: group by (e.g. celltype or genotype)
  #
  # Returns:
  #   p: a Seurat object
  
  DefaultAssay(obj) <- "RNA"
  dat <- FetchData(obj, 
                   vars = c(gene.name, "ident", "genotype", "sample.name"))
  colnames(dat) <- c("Expression", "ident", "genotype", "sample.name")
  
  p <- ggplot(dat, aes_string(x=group.by, y="Expression", col=group.by, 
                              fill=group.by)) + 
    geom_violin(alpha=0.8) +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "point", width=0.5, col="black") +
    xlab("") +
    ylab(paste(gene.name, "expression")) + 
    NoLegend() +
    scale_fill_manual(values=pal)
  p
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., 
                                                 draw_quantiles = NULL) {
                             data <- transform(data, 
                                               xminv = x - violinwidth * (x - xmin), 
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(
                               transform(data, 
                                         x = if (grp %% 2 == 1) xminv else xmaxv), 
                               if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, 
                                              newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & 
                                 !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), 
                                         all(draw_quantiles <= 1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, 
                                                                                    draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), 
                                                  setdiff(names(data), c("x", "y")), 
                                                  drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", 
                                                grid::grobTree(GeomPolygon$draw_panel(newdata, ...), 
                                                               quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", 
                                                GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", 
                              position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, 
                              scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, 
                      draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


GetBgGenes <- function(obj, genes.to.test, dge.results) {
  # Gets a list of background genes expressed in celltypes tested in DGE
  # 
  # Args:
  #   obj: seurat object
  #   genes.to.test: list of genes used for DGE testing 
  #   dge.results: dge results as df
  #
  # Returns:
  #   bg: background genes (a vector)
  
  DefaultAssay(obj) <- "RNA"
  
  cell.idents <- FetchData(obj, vars="ident")
  cell.idents <- subset(cell.idents, ident %in% unique(dge.results$celltype))
  
  bg <- obj[, rownames(cell.idents)]
  bg <- rownames(bg) 
  bg <- bg[bg %in% genes.to.test]
  
  return(bg)
}

CreateContingency <- function(df, gene.list, fdr.cutoff, log2fc.cutoff=Inf) {
  # Create contingency table specifically for DOWNREGULATED genes
  #
  # Args:
  #   df: dge df 
  #   gene.list: list of genes
  #   fdr.cutoff: FDR cutoff
  #   log2fc.cutoff: log fc cutoff
  #
  # Returns:
  #   cont.table: a contigency table
  df$indicator <- ifelse(df$gene %in% gene.list$gene, 1, 0)
  
  diff <- subset(df, FDR < fdr.cutoff & avg_log2FC < log2fc.cutoff)
  nondiff <- subset(df, gene %!in% diff$gene)
  
  # differentially expressed genes in gene list
  diff.in.list <- subset(diff, indicator == 1)
  diff.not.in.list <- subset(diff, indicator == 0)
  
  # non diff
  nondiff.in.list <- subset(nondiff, indicator == 1)
  nondiff.not.in.list <- subset(nondiff, indicator == 0)
  
  cont.table <- matrix(c(nrow(diff.in.list), 
                         nrow(diff.not.in.list), 
                         nrow(nondiff.in.list),
                         nrow(nondiff.not.in.list)), nrow=2)
  return(cont.table)
}


CellTypeEnrichment <- function(dge.df, gene.set, fdr.cutoff=.05, 
                               log2fc.cutoff=0, gene.set.name=NULL) {
  # Tests enirchment of downregulated genes in each celltype; corrects 
  #   p-value via FDR
  celltypes <- unique(dge.df$celltype)
  
  res <- data.frame()
  
  for(i in celltypes) {
    deg <- dge.df %>% filter(celltype == i) %>%
      CreateContingency(gene.set, fdr.cutoff, log2fc.cutoff)
    fet.res <- fisher.test(deg)
    curr <- data.frame(celltype=i, p=fet.res$p.value, OR = fet.res$estimate[[1]])
    res <- rbind(res, curr)
  }
  res$FDR <- p.adjust(res$p, method="BH", n = 9*3)
  
  if(!is.null(gene.set.name)) {
    res$geneset <- gene.set.name
  }
  return(res)
}


PlotEnrichment <- function(dge.df, epi.genes, sfari.genes, dd.genes, 
                           fdr.cutoff = 0.1, log2fc.cutoff = -.1,
                           color.pal) {
  
  epi.res <- CellTypeEnrichment(dge.df, epi.genes, fdr.cutoff, log2fc.cutoff,
                                "Epilepsy genes")
  sfari.res <- CellTypeEnrichment(dge.df, sfari.genes, fdr.cutoff, log2fc.cutoff,
                                  "SFARI genes")
  dd.res <- CellTypeEnrichment(dge.df, dd.genes, fdr.cutoff, log2fc.cutoff,
                               "DD genes")
  
  df <- rbind(epi.res, sfari.res, dd.res)
  df$geneset <- factor(df$geneset, levels=c("DD genes", 
                                            "Epilepsy genes", 
                                            "SFARI genes"))
  
  p <- ggplot(df, aes(x=reorder(celltype, -FDR), 
                      y=-log10(FDR), fill=celltype)) + 
    geom_bar(stat="identity", position="dodge", size=0.25, col="black") + 
    coord_flip() + 
    geom_hline(yintercept = -log10(.05), linetype = "dashed") +
    scale_fill_manual(values = color.pal) +
    facet_grid(cols = vars(geneset)) + 
    NoLegend() +
    xlab("")
  p
}



PlotNumGene <- function(dge.df, epi.genes, sfari.genes, dd.genes, color.pal) {
  # Plot number of epilepsy, sfari, and ddd genes that were tested per celltype
  #   during DGE
  #
  # Args:
  #   dge.df: dge results
  #   epi.genes: df epilepsy genes
  #   sfari.genes: df of sfari genes
  #   dd.genes: df of ddd genes
  #   color.pal: color palette for plotting
  #
  # Returns:
  #   p: a bar chart

  dge.df$epi.gene <- ifelse(dge.df$gene %in% epi.genes$gene, 1, 0)
  dge.df$sfari.gene <- ifelse(dge.df$gene %in% epi.genes$gene, 1, 0)
  dge.df$dd.gene <- ifelse(dge.df$gene %in% dd.genes$gene, 1, 0)
  
  plot.df <- dge.df %>% group_by(celltype) %>% 
    summarize(n.epi = sum(epi.gene), 
              n.sfari = sum(sfari.gene), 
              n.dd = sum(dd.gene))
  plot.df <- melt(plot.df)
  p <- ggplot(plot.df, aes(x=celltype, y=value, fill=celltype)) + 
    geom_bar(stat="identity") +
    facet_grid(cols=vars(variable)) +
    # RotatedAxis() + 
    NoLegend() + 
    scale_fill_manual(values=color.pal) +
    xlab("") +
    ylab("Number of detected genes") +
    coord_flip()
  p
}


PerformGO <- function(dge, logfc.cutoff, fdr.cutoff=.01, celltypes=NULL) {
  # Perform GO analysis via gprofiler
  #
  # Args:
  #   dge: dataframe of dge results
  #   logfc.cutoff: cutoff for sig genes
  #   fdr.cutoff: fdr cutoff for sig genes
  #   celltype: celltype of interest
  #
  # Returns:
  #   res: GO enrichhment results
  if(!is.null(celltypes)) {
    dge <- dge %>% filter(celltype %in% celltypes)
  }
  
  
  if(logfc.cutoff > 0) {
    sig <- dge %>% filter(FDR < fdr.cutoff & avg_log2FC  > logfc.cutoff)
  } else {
    sig <- dge %>% filter(FDR < fdr.cutoff & avg_log2FC  < logfc.cutoff)
  }
  
  res <- gost(sig$gene, correction_method = "fdr", domain_scope = "custom", 
              custom_bg = dge$gene, evcodes = F,
              sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"))
  
  return(res$result)
}


PlotGO <- function(go.res, selected.terms) {
  # Creates barchart of go enrichment results
  #
  # Args:
  #   go.res: go results
  #   selected.terms: terms to include in plot
  #
  # Returns:
  #   p: ggplot barchart
  plot.df <- go.res %>% filter(term_name %in% selected.terms)
  
  p <- ggplot(plot.df, aes(x=reorder(term_name, -p_value), 
                           y=-log10(p_value), fill=-log10(p_value))) + 
    geom_bar(stat="identity") + 
    coord_flip() + 
    xlab("") + 
    geom_hline(yintercept = -log10(.05), linetype = "dashed") + 
    NoLegend() +
    ylab("-log10(FDR)")
  
  p
}


PlotClipTags <- function(clip.data, normalize=F) {
  # Plots number of clip tags per gene; adds line for MEF2C
  # 
  # Args:
  #   clip.data: path to clip seq data
  #
  # Returns:
  #   p: ECDF
  
  clip  <- fread("data/clip_tag_tallies.txt")
  clip <- clip %>% filter(n_clip_tags > 1)
  
  x.label <- "Number of Hnrnpu binding sites per gene"
  
  if(normalize) {
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    x <- getBM(attributes = c("mgi_symbol", "transcript_length"),  
               mart=mouse)
    y <- x %>% group_by(mgi_symbol) %>% top_n(1, transcript_length)
    
    clip <- left_join(clip, y, by=c('gene' = 'mgi_symbol'))
    
    clip$n_clip_tags <- clip$n_clip_tags / clip$transcript_length
    x.label <- "Number of Hnrnpu binding sites per base per gene"
  }
  
  mef2c <- clip %>% filter(gene == "Mef2c")
  mef2c.n <- mef2c$n_clip_tags
  
  p <- clip %>% filter((gene %in% rownames(hippo) | 
                          gene %in% rownames(ctx))) %>% 
    ggplot(aes(x=n_clip_tags)) +
    stat_ecdf(col = "steelblue") + 
    geom_vline(xintercept = mef2c.n, linetype = "dashed") + 
    ylab("Cumulative proportion of genes") + 
    xlab(x.label)
  
  p
}


PlotCMAPCandidates <- function(df, plot.title=NULL, n=20) {
  # Plots top n CMAP compounds
  #
  # Args:
  #   df: df containing CMAP results
  #   plot.title: plot title
  #   n: number of compounds to plot
  #
  # Returns:
  #   p: dot plot
  
  pal <- colorRampPalette(rev(brewer.pal(8, "YlGnBu")))(n)
  
  p <- df %>% 
    top_n(Score, n = -n) %>% 
    ggplot(aes(x=reorder(Name, -Score), y=Score, fill=Score)) + 
    geom_point(stat="identity", shape = 21, stroke=0.2, size=2) + 
    ylim(-97, -100) +
    coord_flip() + 
    scale_fill_gradientn(colours = pal) + 
    NoLegend() + 
    xlab("Compound") + 
    ylab("Connectivity Score") + 
    ggtitle(plot.title) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p
}


PlotCMAPDist <- function(df, plot.title=NULL) {
  # Plot histogram of CMAP scores
  #
  # Args:
  #   df: dataframe of CMAP results
  #   plot.title: title of plot
  #
  # Returns:
  #   p: a histogram
  p <- df %>% filter(Score <= 0) %>% 
    ggplot(aes(x=Score)) + 
    geom_histogram(col="black", size = 0.1, binwidth = 4, fill="lightblue") +
    NoLegend() +
    xlim(0, -100) + 
    xlab("Connectivity Score") + 
    ylab("Count") + 
    ggtitle(plot.title) + 
    geom_vline(xintercept = -90, linetype = "dashed") +
    theme(plot.title = element_text(hjust = 0.5))
  
  p
}

################################################################################
# Main                                                                         #
################################################################################
# Figure 1: clusters -----------------------------------------------------------
# Umap plots
hippo.umap <- VisualizeUMAP(hippo) 
ctx.umap <- VisualizeUMAP(ctx)

# Heatmap / dot plots for marker viz
hippo.dotplot <- CreateDotPlot(hippo, "data/hippocampal_markers.txt")
ctx.dotplot <- CreateDotPlot(ctx, "data/cortical_markers.txt")

# Hnrnpu expression
hnrnpu.hippo <- hippo %>% subset(genotype == "WT") %>% 
  JitterPlot("Hnrnpu", "ident") + 
  scale_color_manual(values = pal) +
  # RotatedAxis() + 
  coord_flip() + 
  ylim(0, 3)

hnrnpu.ctx <- ctx %>% subset(genotype == "WT") %>% 
  JitterPlot("Hnrnpu", "ident") + 
  scale_color_manual(values = pal) +
  # RotatedAxis() +
  coord_flip() + 
  ylim(0, 3)

fig1 <- plot_grid(hippo.umap, hippo.dotplot,  hnrnpu.hippo,
                  ctx.umap, ctx.dotplot,
                  hnrnpu.ctx, 
                  ncol = 3, 
                  rel_widths = c(1, 1, 0.65),
                  labels="AUTO")

# View number of cells per cluster
PlotCellCounts(hippo, "genotype") + geom_hline(yintercept=300)
PlotCellCounts(ctx, "genotype") + geom_hline(yintercept=300)

save_plot("figures/fig1.pdf", fig1, base_width = 11, base_height = 6.5)


# Figure 2: DGE ----------------------------------------------------------------
# Differential gene expression testing
DefaultAssay(hippo) <- "RNA"
DefaultAssay(ctx) <- "RNA"

# Uncomment to run DGE..........................................................
# 
# genes.to.test <- fread("data/genes_for_dge.txt", header=F)
# genes.to.test <- subset(genes.to.test, genes.to.test$V1 %in% rownames(hippo))
# genes.to.test <- genes.to.test$V1
# 
# hippo.dge <- MultiDGE(hippo, min.pct = 0.1, logfc.threshold = 0,
#                       out.file = "results/dge/hippocampus_dge.txt",
#                       genes = genes.to.test,
#                       exclude.celltypes = c("Prolif.", "Astrocyte",
#                                             "Cajal", "OPC"))
# 
# genes.to.test <- fread("data/genes_for_dge.txt", header=F)
# genes.to.test <- subset(genes.to.test, genes.to.test$V1 %in% rownames(hippo))
# genes.to.test <- genes.to.test$V1

# ctx.dge <- MultiDGE(ctx, min.pct = 0.1, logfc.threshold = 0, 
#                     out.file = "results/dge/cortex_dge.txt", 
#                     genes = genes.to.test,
#                     exclude.celltypes = c("Prolif.", "Astrocyte", 
#                                           "Cajal", "OPC"))

# Downsampled DGE
# hippo.ds.dge <- MultiDGE(hippo, min.pct = 0.1, logfc.threshold = 0, 
#                          out.file = "results/dge/hippocampus_dge_downsampled.txt", 
#                          genes = genes.to.test, 
#                          max.cells = 300,
#                          exclude.celltypes = c("Prolif.", "Astrocyte", 
#                                                "Cajal", "OPC", "Radial Glia"), 
#                          test.use = "MAST_downsample")
# 
# ctx.ds.dge <- MultiDGE(ctx, min.pct = 0.1, logfc.threshold = 0, 
#                        out.file = "results/dge/cortex_dge_downsampled.txt", 
#                        genes = genes.to.test, 
#                        max.cells = 300,
#                        exclude.celltypes = c("Prolif.", "Astrocyte", "Cajal", 
#                                              "OPC", "Radial Glia", "LGE"), 
#                        test.use = "MAST_downsample")

hippo.dge <- fread("results/dge/hippocampus_dge.txt")
hippo.ds.dge <- fread("results/dge/hippocampus_dge_downsampled.txt")
hippo.dge$FDR <- p.adjust(hippo.dge$p_val, method = "BH")
hippo.ds.dge$FDR <- p.adjust(hippo.ds.dge$p_val, method = "BH")

ctx.dge <- fread("results/dge/cortex_dge.txt")
ctx.ds.dge <- fread("results/dge/cortex_dge_downsampled.txt")
ctx.dge$FDR <- p.adjust(ctx.dge$p_val, method = "BH")
ctx.ds.dge$FDR <- p.adjust(ctx.ds.dge$p_val, method = "BH")

hippo.volcano <- VolcanoPlot(hippo.dge, log2fc.cutoff = .14, 
                             label.log2fc = .4, label.fdr = 10^-10, 
                             color.pal = pal) +
  xlim(c(-1.3, 1.3)) + 
  ylim(c(0, 50))

ctx.volcano <- VolcanoPlot(ctx.dge, log2fc.cutoff = .14, 
                           label.log2fc = .35, label.fdr = 10^-8, 
                           color.pal = pal) + 
  ylim(c(0, 50)) + 
  xlim(c(-1.3, 1.3))


hippo.burden <- PlotNumDEG(hippo.ds.dge, color.pal = pal, 
                           log2fc.cutoff = .14, fdr.cutoff = .05) +
  ylim(0, 35) + 
  coord_flip()

ctx.burden <- PlotNumDEG(ctx.ds.dge, color.pal = pal, 
                         log2fc.cutoff = .14, fdr.cutoff = .05) + 
  ylim(0,35) + 
  coord_flip()

# GO enrichment
hippo.go <- PerformGO(hippo.dge, logfc.cutoff = -.14, fdr.cutoff = .05)
hippo.go.up <- PerformGO(hippo.dge, logfc.cutoff = .14, fdr.cutoff = .05)
fwrite(hippo.go, "results/GO/hippo_down.csv")
fwrite(hippo.go.up, "results/GO/hippo_up.csv")

hippo.terms <- c("neuron projection development", "neuron differentiation", 
                 "axonogenesis", "movement of cell or subcellular component",
                 "locomotion", "neuron migration", 
                 "positive regulation of gene expression")
hippo.go.plot <-PlotGO(hippo.go, hippo.terms) +
  ylim(0, 10)


ctx.go <- PerformGO(ctx.dge, logfc.cutoff = -.14, fdr.cutoff = .05)
ctx.go.up <- PerformGO(ctx.dge, logfc.cutoff = .14, fdr.cutoff = .05)
fwrite(ctx.go, "results/GO/ctx_down.csv")
fwrite(ctx.go.up, "results/GO/ctx_up.csv")

ctx.terms <- c("regulation of transcription by RNA polymerase II", 
               "regulation of RNA metabolic process", 
               "neuron differentiation",
               "developmental process",
               "axon development",
               "neurogenesis",
               "synaptic transmission, glutamatergic")
ctx.go.plot <-PlotGO(ctx.go, ctx.terms) +
  ylim(0, 10)

# Gene list enrichments
# enrichment
epi.genes <- fread("data/gene_lists/epi25_mouse.txt")
sfari.genes <- fread("data/gene_lists/sfari_mouse.txt")
dd.genes <- fread("data/gene_lists/ddd_brain_mouse.txt")

hippo.enrich <- PlotEnrichment(hippo.dge, epi.genes, sfari.genes, dd.genes, 
                               log2fc.cutoff = -.074, fdr.cutoff = 0.1,
                               color.pal = pal)

ctx.enrich <- PlotEnrichment(ctx.dge, epi.genes, sfari.genes, dd.genes, 
                             log2fc.cutoff = -.074, fdr.cutoff = 0.1, 
                             color.pal = pal)


enrich.plot <- plot_grid(hippo.enrich, ctx.enrich, nrow=2)

row1 <- plot_grid(hippo.volcano, hippo.burden,  
          ctx.volcano, ctx.burden, 
          rel_widths = c(1, 0.7, 1, 0.7),
          nrow=1, labels = "AUTO",
          label_size = 10)

row2 <- plot_grid(hippo.go.plot, ctx.go.plot, enrich.plot, 
                  nrow=1, rel_widths = c(0.8, 0.8, 1), 
                  labels = c("E", "F", "G"), 
                  label_size = 10)

fig2 <- plot_grid(row1, row2, nrow=2)  
save_plot("figures/fig2.pdf", fig2, base_width = 11, base_height = 6.5)

# Figure 3: MEF2C --------------------------------------------------------------
hippo.mef2c <- SplitVlnPlot(hippo, "Mef2c", 
                            exclude.celltypes = c("Prolif.", "Astrocyte", 
                                                  "Cajal", "OPC")) + 
  ggtitle("Mef2c expression (hippocampus)") + 
  theme(legend.position = c(0.02, 1), legend.justification = c(0.02, 1))


ctx.mef2c <- SplitVlnPlot(ctx, "Mef2c", 
                          exclude.celltypes = c("Prolif.", "Astrocyte", 
                                                "Cajal", "OPC")) + 
  ggtitle("Mef2c expression (cortex)") +
  theme(legend.position = c(0.02, 1), legend.justification = c(0.02, 1))


clip.cdf <- PlotClipTags("data/clip_tag_tallies.txt")
clip.normalized.cdf <- PlotClipTags("data/clip_tag_tallies.txt", normalize = T)

fig3 <- plot_grid(hippo.mef2c, ctx.mef2c, clip.cdf, clip.normalized.cdf, nrow=2,
                  labels="AUTO", label_size = 10, 
                  rel_heights = c(1, 0.9))

save_plot("figures/fig3.pdf", fig3, base_width = 6.5, base_height = 5)


# Figure4: CMAP analysis -------------------------------------------------------
# Run pseudo-bulk analysis
hippo.bulk.dge <- BulkDGE(hippo, genes.to.test)
write.table(hippo.bulk.dge, "results/dge/hippocampus_pseudo_bulk.txt", 
            quote=F, sep="\t", row.names = F)

subiculum.cmap <- fread("results/cmap/subiculum_downregulated.txt")
bulk.cmap <- fread("results/cmap/bulk_downregulated.txt")

p1 <- PlotCMAPDist(subiculum.cmap, 
                   plot.title = "Subiculum signature reversal")

p2 <- PlotCMAPDist(bulk.cmap, 
                   plot.title = "Pseudo-bulk hippocampus signature reversal")

p3 <- PlotCMAPCandidates(subiculum.cmap, plot.title = "Top subiculum candidates")
p4 <- PlotCMAPCandidates(bulk.cmap, plot.title = "Top pseudo-bulk candidates")

fig4 <- plot_grid(p1, p2, p3, p4, nrow=2, labels = "AUTO", 
                  label_size = 10)


save_plot("figures/fig4.pdf", fig4, base_width = 6.5, base_height = 6)


################################################################################
# Supplementary figures                                                        #
################################################################################
hippo.umap.split <- hippo.umap + 
  facet_wrap(~sample.name)

ctx.umap.split <- VisualizeUMAP(ctx) + 
  facet_wrap(~sample.name)

figS4 <- plot_grid(hippo.umap.split, ctx.umap.split, nrow=2, labels = "AUTO", label_size = 8)
save_plot("figures/figureS4.pdf", figS4, base_height = 8, base_width = 5)

# num cells
n.hippo <- PlotCellCounts(hippo)
n.ctx <- PlotCellCounts(ctx)

figS5 <- plot_grid(n.hippo, n.ctx, nrow = 1, labels = "AUTO")
save_plot("figures/figureS5.pdf", figS5, base_height = 4, base_width = 8)



PlotNumGene(hippo.dge, epi.genes, sfari.genes, dd.genes, color.pal = hippo.pal)
PlotNumGene(ctx.dge, epi.genes, sfari.genes, dd.genes, color.pal = ctx.pal)

