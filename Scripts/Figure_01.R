#Figure 1 Components
#Final tweaks and labels were applied in Adobe Illustrator
#####

# Day Samples TSNE Plot

day_tsneplot <- TSNEPlot(Pineal2_Day1n2, do.return = T, pt.size = 0.10) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                                                     axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.border = element_blank())
# Violin Plot Grid 

object <- Pineal2_Day1n2
celltype.labels <- c("aPinealocyte","bPinealocyte","aAstrocyte","bAstrocyte","gAstrocyte","aMicroglia","bMicroglia", "VLMCs","Endothelial")
violin.palette <- c("#E66AA8","#B37BB5","#F3766E","#D29329","#92AB3C","#2AB34B","#2EBA9C","#0BB9E3","#7094CD")
genes.use <- c("Asmt","Cngb1","Penk","S100b","Lum","Dcn","Lyz2","C1qa","RT1-Bb","Vwf","Emcn")

celltypes <- object@ident
num.cells <- dim(object@data)[2]
violin.matrix <- as.data.frame(matrix(nrow=num.cells*length(genes.use), ncol=3))
colnames(violin.matrix) <- c("expression","celltype","gene")
violin.matrix$celltype <- celltypes

for(i in 1:length(genes.use)){
  violin.matrix$expression[((i*num.cells)-num.cells+1):(i*num.cells)] <- object@data[genes.use[i],]
  violin.matrix$gene[((i*num.cells)-num.cells):(i*num.cells)] <- genes.use[i]
}

violin.matrix$gene <- factor(violin.matrix$gene, levels = genes.use)
violin.matrix$celltype <- factor(violin.matrix$celltype, levels(violin.matrix$celltype)[c(9,8,1,2,3,4,5,6,7)])

violin.grid <- ggplot(data = violin.matrix, aes(y=expression, x="")) + geom_violin(aes(fill=celltype)) + facet_grid(gene ~ celltype, switch = 'both') + 
  ylab(NULL) + xlab(NULL) + scale_y_discrete(labels = celltype.labels) +
  theme(panel.spacing=unit(0.01, "lines"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        legend.position = 'none', axis.ticks.y = element_line(size=0), 
        axis.title.x = element_text(size=10), axis.text.x = element_text(size = 10), 
        strip.text.y = element_blank(), 
        strip.text.x = element_blank(), 
        strip.background = element_rect(fill = 'white', size=0)) + scale_fill_manual(values=violin.palette)

# Dendrogram 

dendrogram <- BuildClusterTree(Pineal2_Day1n2)
dendrogram <- dendrogram@cluster.tree[[1]]
dendrogram$tip.label <- c("aPinealocyte","bPinealocyte","aAstrocyte","bAstrocyte","gAstrocyte","aMicroglia","bMicroglia", "VLMCs","Endothelial")
#dendrogram$tip.label <- c("","","","","","","","","")

# Output

svglite::svglite(file = "./day_tsne_plot.svg", width = 3.5, height = 3.5)
day_tsneplot
dev.off()

pdf(file = "./violin_grid.pdf", width = 3.2, height = 4.5)
violin.grid
dev.off()

png(filename = "./day_dendrogram.png", width = 3.5, height = 5, res = 300, units = "in")
ape::plot.phylo(dendrogram)
dev.off()