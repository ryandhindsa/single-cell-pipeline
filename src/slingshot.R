library(slingshot)
library(RColorBrewer)


obj <- readRDS("~/Dropbox/scRNA/projects/hnrnpu_cortex/results/hnrnpu_cortex_integrated.rds")
wt <- subset(obj, genotype %in% c("WT"))

sce <- as.SingleCellExperiment(wt, assay = "RNA")
sce@reducedDims$PCA2 <- sce@reducedDims$PCA[, 1:2]
sce <- slingshot(sce, clusterLabels = "celltype", start.clus = "NSC", reducedDim = "PCA2")

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce)$PCA2, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

plot(reducedDims(sce)$PCA2, col = brewer.pal(9,'Set1')[sce$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', show.constraints = TRUE)
FeaturePlot(obj, reduction = "pca", "GAD2", min.cutoff = "q9", dims = 1:2)

###############################################################################
