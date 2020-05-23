# Imports ----------------------------------------------------------------------
library(ggplot2)
library(cowplot)
library(tidyverse)

# Functions --------------------------------------------------------------------
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
  # gene.set <- filter(gene.set, gene != "Hnrnpu")
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

# Main -------------------------------------------------------------------------
hippo.dge <- fread("results/dge/hippocampus_dge.txt")
ctx.dge <- fread("results/dge/cortex_dge.txt")

# enrichment
epi.genes <- fread("data/gene_lists/epi25_mouse.txt")
sfari.genes <- fread("data/gene_lists/sfari_mouse.txt")
dd.genes <- fread("data/gene_lists/ddd_brain_moust.txt")

PlotEnrichment(hippo.dge, epi.genes, sfari.genes, dd.genes, 
               log2fc.cutoff = -.074, fdr.cutoff = 0.1,
               color.pal = hippo.pal)

PlotEnrichment(ctx.dge, epi.genes, sfari.genes, dd.genes, 
               log2fc.cutoff = -.074, fdr.cutoff = 0.1, 
               color.pal = ctx.pal)
