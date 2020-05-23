props <- prop.table(x = table(Idents(obj), obj@meta.data[, "sample.name"]), 
                    margin = 2)
x <- data.frame(props)

ggplot(x, aes(x=Var1,  y=Freq, fill=Var2)) + 
  geom_bar(stat="identity", position="dodge", col="black") + 
  RotatedAxis()

table(Idents(obj), obj@meta.data[, "sample.name"])

obj <- readRDS("~/Dropbox/scRNA/projects/map1b_flna_co53/results/map1b_flna_co53_integrated.rds")
DimPlot(obj)
unique(obj$sample.name)
DimPlot(obj, label=T)
FeaturePlot(obj, "EOMES", min.cutoff = "q9")
DimPlot(obj)
obj <- RenameIdents(obj, "1" = "Neuron")
FeaturePlot(obj, "TCF7L2", min.cutoff = "q9")
DimPlot(obj, split.by = "sample.name")
