#Figure 4 Components
#Final tweaks and labels were applied in Adobe Illustrator
#####

# Rename Celltypes
levels(Pineal2_Day1n2@ident) <- c("aPinealocyte","bPinealocyte","aAstrocyte","bAstrocyte","gAstrocyte","aMicroglia","bMicroglia", "VLMCs","Endothelial")
levels(Pineal2_Combined.2@ident) <- c("aPinealocyte","bPinealocyte","aAstrocyte","bAstrocyte","gAstrocyte","aMicroglia","bMicroglia", "VLMCs","Endothelial")

# Comparing Astrocyte Subgroups
#Find markers using Wilcoxon Ranked Sum Test
object <- Pineal2_Day1n2
ident.1 <- "aAstrocyte"
ident.2 <- "bAstrocyte"
ident.3 <- "gAstrocyte"

marker.report.raw.1 <- FindMarkers(object, ident.1 = ident.1, ident.2 = c(ident.2,ident.3), logfc.threshold = 0.25, min.pct = 0.10, test.use = 'wilcox', pseudocount.use = 0.01)
marker.report.raw.2 <- FindMarkers(object, ident.1 = ident.2, ident.2 = c(ident.1,ident.3), logfc.threshold = 0.25, min.pct = 0.10, test.use = 'wilcox', pseudocount.use = 0.01)
marker.report.raw.3 <- FindMarkers(object, ident.1 = ident.3, ident.2 = c(ident.1,ident.2), logfc.threshold = 0.25, min.pct = 0.10, test.use = 'wilcox', pseudocount.use = 0.01)
marker.report.astro.a <- GenerateMarkerReport(marker.report.raw = marker.report.raw.1, object = object, ident.1 = ident.1, ident.2 = c(ident.2,ident.3))
marker.report.astro.b <- GenerateMarkerReport(marker.report.raw.2, object = object, ident.1 = ident.2, ident.2 = c(ident.1,ident.3))
marker.report.astro.c <- GenerateMarkerReport(marker.report.raw.3, object = object, ident.1 = ident.3, ident.2 = c(ident.1,ident.2))
marker.list <- list()

marker.list$ident.1 <- (marker.report.astro.a %>% filter(p.val.fdr<0.05, ident.2.up==F, ident.1.fc>=2.0) %>% arrange(desc(effect.size)))$gene
marker.list$ident.2 <- (marker.report.astro.b %>% filter(p.val.fdr<0.05,ident.2.up==F, ident.1.fc>=2.0) %>% arrange(desc(effect.size)))$gene
marker.list$ident.3 <- (marker.report.astro.c %>% filter(p.val.fdr<0.05,ident.2.up==F, ident.1.fc>=2.0) %>% arrange(desc(effect.size)))$gene
marker.list$combined <- c(marker.list$ident.1[1:6], marker.list$ident.2[1:6], marker.list$ident.3[c(1:6)])

object.sub <- SubsetData(object, ident.use = c("aAstrocyte","bAstrocyte","gAstrocyte"))
object.sub <- ScaleData(object.sub, display.progress = F)
cells.use <- c(sample(x = WhichCells(object.sub, ident = "aAstrocyte"), size = 100), WhichCells(object.sub, ident = c("bAstrocyte","gAstrocyte")))

#Top 6 markers of each cell type, ranked by effect size
DoHeatmap(object.sub, use.scaled = T, genes.use = marker.list$combined, cells.use = cells.use, slim.col.label = T, col.low = "#3066BE", col.high = "#D81139", remove.key = F,group.label.rot = F, group.cex = 10, rotate.key = T)

#Output
pdf(file = "./AstrocyteHeatmap_nokey_small.pdf", width = 3.54, height = 2.36)
DoHeatmap(object.sub, use.scaled = T, genes.use = marker.list$combined, cells.use = cells.use, slim.col.label = T, col.low = "#3066BE", col.high = "#D81139", remove.key = T)
dev.off()