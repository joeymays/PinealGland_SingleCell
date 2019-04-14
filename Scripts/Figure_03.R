#Figure 3 Components
#Final tweaks and labels were applied in Adobe Illustrator
#####

#Rename Celltypes 
levels(Pineal2_Day1n2@ident) <- c("aPinealocyte","bPinealocyte","aAstrocyte","bAstrocyte","gAstrocyte","aMicroglia","bMicroglia", "VLMCs","Endothelial")
levels(Pineal2_Combined.2@ident) <- c("aPinealocyte","bPinealocyte","aAstrocyte","bAstrocyte","gAstrocyte","aMicroglia","bMicroglia", "VLMCs","Endothelial")

#Comparing Subgroups
#Find markers using Wilcoxon Ranked Sum Test
object <- Pineal2_Day1n2
ident.1 <- "aPinealocyte"
ident.2 <- "bPinealocyte"

marker.report.raw <- FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0.25, min.pct = 0.15, test.use = 'wilcox', pseudocount.use = 0.01)


#Generate Marker Report, including FDR at 0.05, Fold Change, Effect Size
genes.use <- rownames(marker.report.raw)
cells.ident.2 <- WhichCells(object, ident.2)
cells.ident.1 <- WhichCells(object, ident.1)
cells.ident.2.n <- length(cells.ident.2)
cells.ident.1.n <- length(cells.ident.1)

mean.ident.2 <- apply(X = object@data[genes.use, cells.ident.2, drop = F], 
                      MARGIN = 1, FUN = function(x) mean(x = expm1(x = x)))

mean.ident.1 <- apply(X = object@data[genes.use, cells.ident.1, drop = F], 
                      MARGIN = 1, FUN = function(x) mean(x = expm1(x = x)))

sd.ident.2 <- apply(X = object@data[genes.use, cells.ident.2, drop = F], 
                    MARGIN = 1, FUN = function(x) sd(x = expm1(x = x)))

sd.ident.1 <- apply(X = object@data[genes.use, cells.ident.1, drop = F], 
                    MARGIN = 1, FUN = function(x) sd(x = expm1(x = x)))

pooled.sd <- sqrt((((cells.ident.1.n-1)*(sd.ident.1^2))+((cells.ident.2.n-1)*(sd.ident.2^2)))/(cells.ident.1.n + cells.ident.2.n - 2))

marker.report <- as.data.frame(matrix(nrow = length(rownames(marker.report.raw)), ncol = 0, dimnames = list(genes.use,c())))

marker.report$pct.1 <- marker.report.raw[,3]
marker.report$pct.2 <- marker.report.raw[,3]
marker.report$p.val.fdr <- p.adjust(marker.report.raw[,1], method = "fdr")
marker.report$effect.size <- abs((mean.ident.2 - mean.ident.1)/(pooled.sd))
marker.report$mean.1 <- mean.ident.1
marker.report$mean.2 <- mean.ident.2
marker.report$log.fc <- marker.report.raw[,2]
marker.report$ident.1.fc <- (marker.report$mean.1+0.01)/(marker.report$mean.2+0.01)
marker.report$ident.2.fc <- (marker.report$mean.2+0.01)/(marker.report$mean.1+0.01)
marker.report$ident.2.up <- (marker.report.raw[,2])<0
marker.report$gene <- rownames(marker.report.raw)


#Heatmap
#Select 10 highest genes by effect size; Samples of 250 cells
marker.list <- list()

marker.list$ident.1 <- (marker.report %>% filter(p.val.fdr<0.05, ident.2.up==F, ident.1.fc>=2.0) %>% arrange(desc(effect.size)))$gene
marker.list$ident.2 <- (marker.report %>% filter(p.val.fdr<0.05,ident.2.up==T, ident.2.fc>=2.0) %>% arrange(desc(effect.size)))$gene
marker.list$combined <- c(marker.list$ident.1[1:10], c(marker.list$ident.2[c(1:5,7,9,11,12,13)]))

tph1.sample <- WhichCells(object, "aPinealocyte")
tph1.sample <- sample(x = tph1.sample, size = 250, replace = F)
asmt.sample <- WhichCells(object, "bPinealocyte")
asmt.sample <- sample(x = asmt.sample, size = 250, replace = F)

cells.use <- c(asmt.sample, tph1.sample)

object.sub <- SubsetData(object, ident.use = c("bPinealocyte", "aPinealocyte"))
object.sub <- ScaleData(object.sub, display.progress = F)


#Top 10 markers of each cell type, ranked by effect size
pinealocyte_hmap <- DoHeatmap(object.sub, use.scaled = T, cells.use = cells.use, genes.use = marker.list$combined, 
               slim.col.label = T, col.low = "#3066BE", col.high = "#D81139", remove.key = F)

# Violin Plots

#All DE mitochondrial genes (p < 0.05, fc >= 2.0)
mito.genes <- (marker.report[grep("Mt-", marker.report$gene),] %>% filter(p.val.fdr < 0.05, ident.1.fc >= 2.0) %>% arrange(desc(effect.size)))$gene
#Top 20 DE ribosomal genes by effect size (p < 0.05, fc >= 2.0)
ribo.genes <- (marker.report[grep("Rp[sl]", marker.report$gene),] %>% filter(p.val.fdr < 0.05, ident.2.fc >= 2.0) %>% arrange(desc(effect.size)))$gene[1:20]
asmt.genes <- c("Asmt")
ggamma.genes <- c("Gngt1","Gngt2","Gng10","Gng13")
genes.use.list <- list(mito.genes, ribo.genes, asmt.genes, ggamma.genes)
names(genes.use.list) <- c("mito.genes", "ribo.genes", "asmt.genes", "ggamma.genes")

celltypes.use <- c("aPinealocyte", "bPinealocyte")
object <- Pineal2_Day1n2
violin.data <- list()
cells.use.1 <- WhichCells(object, ident = celltypes.use[1])
cells.use.2 <- WhichCells(object, ident = celltypes.use[2])

for(i in 1:length(genes.use.list)){
  to.plot <- t(as.matrix(expm1(object@data[genes.use.list[[i]],cells.use.1])))
  to.plot <- reshape2::melt(to.plot)
  
  to.plot2 <- as.data.frame(matrix(ncol = 0, nrow = length(to.plot$value)))
  to.plot2$value <- to.plot$value
  to.plot2$class <- names(genes.use.list)[i]
  to.plot2$celltype <- celltypes.use[1]
  
  to.plot <- t(as.matrix(expm1(object@data[genes.use.list[[i]],cells.use.2])))
  to.plot <- reshape2::melt(to.plot)
  
  to.plot3 <- as.data.frame(matrix(ncol = 0, nrow = length(to.plot$value)))
  to.plot3$value <- to.plot$value
  to.plot3$class <- names(genes.use.list)[i]
  to.plot3$celltype <- celltypes.use[2]
  
  to.plot.f <- rbind(to.plot2,to.plot3)
  to.plot.f$celltype <- factor(to.plot.f$celltype, levels = rev(c("aPinealocyte", "bPinealocyte")))
  violin.data[[names(genes.use.list)[i]]]<- to.plot.f
}

plot.list <- list()
plot.list$a <- ggplot(violin.data$mito, aes(x=celltype,y=value,fill=celltype)) + geom_violin(scale='width', trim=T) + 
  coord_flip() + stat_summary(fun.y = "mean", geom = "point", shape = 16, size = 3, color = "midnightblue") + 
  theme(legend.position = 'none') + xlab("") + ylab("Counts") + scale_x_discrete(labels = rev(c("aPinealocyte", "bPinealocyte"))) + ggtitle("mito")

plot.list$b <- ggplot(violin.data$ribo, aes(x=celltype,y=log1p(value),fill=celltype)) + geom_violin(scale='width', trim=T) + 
  coord_flip() + stat_summary(fun.y = "mean", geom = "point", shape = 16, size = 3, color = "midnightblue") + 
  theme(legend.position = 'none') + xlab("") + ylab("Log n Counts") + scale_x_discrete(labels = rev(c("aPinealocyte", "bPinealocyte"))) + ggtitle("ribo")

plot.list$c <- ggplot(violin.data$asmt, aes(x=celltype,y=value,fill=celltype)) + geom_violin(scale='width', trim=T) + 
  coord_flip() + stat_summary(fun.y = "mean", geom = "point", shape = 16, size = 3, color = "midnightblue") + 
  theme(legend.position = 'none') + xlab("") + ylab("Counts") + scale_x_discrete(labels = rev(c("aPinealocyte", "bPinealocyte"))) +ggtitle("asmt")

plot.list$d <- ggplot(violin.data$ggamma, aes(x=celltype,y=log1p(value),fill=celltype)) + geom_violin(scale='width', trim=T) + 
  coord_flip() + stat_summary(fun.y = "mean", geom = "point", shape = 16, size = 3, color = "midnightblue") + 
  theme(legend.position = 'none') + xlab("") + ylab("Log n Counts") + scale_x_discrete(labels = rev(c("aPinealocyte", "bPinealocyte"))) + ggtitle("ggamma")

cowplot::plot_grid(plotlist = plot.list, ncol = 2, labels = 'AUTO')

#Summary Stats
for(i in 1:length(names(genes.use.list))){
  cat(paste("\n",names(genes.use.list)[i]))
  print(knitr::kable(violin.data[[i]] %>% group_by(celltype) %>% summarize(mean = mean(value), median = median(value))))
}

#Output
pdf(file = "./PinealocyteHeatmap_nokey_small.pdf", width = 3.54, height = 2.36)
DoHeatmap(object.sub, use.scaled = T, cells.use = cells.use, genes.use = marker.list$combined, 
          slim.col.label = T, col.low = "#3066BE", col.high = "#D81139", remove.key = T)
dev.off()

pdf(file = "./PinealocyteViolins.pdf", width = 8, height = 4)
cowplot::plot_grid(plotlist = plot.list, ncol = 2, labels = 'AUTO')
dev.off()