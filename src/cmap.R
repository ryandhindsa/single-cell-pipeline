sub <- fread("results/cmap/subiculum_downregulated.txt")
bulk <- fread("results/cmap/bulk_downregulated.txt")


sub$Score <- sub$Score * -1
bulk$Score <- bulk$Score * -1

pal <- colorRampPalette(rev(brewer.pal(8, "YlGnBu")))(20)

sub %>% top_n(Score, n = -20) %>% 
  ggplot(aes(x=reorder(Name, -Score), y=Score, fill=Score)) + 
  geom_point(stat="identity", shape = 21, stroke=0.2, size=3) + 
  ylim(-97, -100) +
  coord_flip() + 
  scale_fill_gradientn(colours = pal) + 
  NoLegend() + 
  xlab("Compound") + 
  ylab("Connectivity Score") + 
  ggtitle("Subiculum signature")
  # scale_y_reverse()
  

bulk %>% top_n(Score, n = -20) %>% 
  ggplot(aes(x=reorder(Name, -Score), y=Score, fill=Score)) + 
  geom_point(stat="identity", shape = 21, stroke=0.2, size=3) + 
  ylim(-97, -100) +
  coord_flip() + 
  scale_fill_gradientn(colours = pal) + 
  NoLegend() + 
  xlab("Compound") + 
  ylab("Connectivity Score") + 
  ggtitle("Pseudo-bulk hippocampus signature")


sub.top <- subiculum.cmap %>% filter(Score < - 90)
bulk.top <- bulk.cmap %>% filter(Score < -90)

nrow(subset(sub.top, sub.top$Name %in% bulk.top$Name))
nrow(subset(sub.top, sub.top$Name %in% bulk.top$Name)) / nrow(sub.top)

pal <- colorRampPalette((brewer.pal(8, "Blues")))(100)

bulk.cmap$candidate <- ifelse(subiculum.cmap$Score < -90, "T", "F")
bulk.cmap %>% filter(Score <= 0) %>% 
  ggplot(aes(x=Score)) + 
  geom_histogram(col="black", size = 0.1, binwidth = 4, fill="lightblue") +
  geom_vline(xintercept = -90, linetype = "dashed") +
  NoLegend() +
  xlim(0, -100)
  # scale_fill_manual(values = pal)

