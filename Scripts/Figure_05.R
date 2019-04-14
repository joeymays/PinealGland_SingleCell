#Figure 5 Components
#Final tweaks and labels were applied in Adobe Illustrator
#####

#Generate Differential Expression Report
object <- Pineal2_Combined.2
celltypes <- levels(object@ident)
test.use <- "wilcox" 
de.test.genes <- rownames(Pineal2_Combined.2@data)

for(i in celltypes){
  object <- SetIdent(object, cells.use = WhichCells(object, ident = i, subset.name = "treatment", accept.value = "Day"), ident.use = paste(i,".day",sep=""))
  object <- SetIdent(object, cells.use = WhichCells(object, ident = i, subset.name = "treatment", accept.value = "Night"), ident.use = paste(i,".night",sep=""))
}

object@ident <- factor(object@ident,levels(object@ident)[c(3,4,15,16,7,8,9,10,11,12,1,2,13,14,5,6,17,18)])

test.celltypes <- c()
for(i in celltypes){
  test.celltypes <- c(test.celltypes,paste(i,".night",sep=""),paste(i,".day",sep=""))
}

celltype.count = 0
night.de.reports <- list()

for(i in celltypes){
  de.test.result <- FindMarkers(object = object, ident.1 = test.celltypes[celltype.count+2], ident.2 = test.celltypes[celltype.count+1], genes.use = de.test.genes, 
                                logfc.threshold = 0, test.use = test.use, min.pct = 0, min.cells.group = 0, pseudocount.use = 0.01)
  
  genes.use = rownames(de.test.result)
  
  cells.treated = WhichCells(object, paste(i,".night",sep=""))
  cells.control = WhichCells(object, paste(i,".day",sep=""))
  cells.treated.n = length(cells.treated)
  cells.control.n = length(cells.control)
  
  mean.treated <- apply(X = object@data[genes.use, cells.treated, drop = F], 
                        MARGIN = 1, FUN = function(x) mean(x = expm1(x = x)))
  
  mean.control <- apply(X = object@data[genes.use, cells.control, drop = F], 
                        MARGIN = 1, FUN = function(x) mean(x = expm1(x = x)))
  
  sd.treated <- apply(X = object@data[genes.use, cells.treated, drop = F], 
                      MARGIN = 1, FUN = function(x) sd(x = expm1(x = x)))
  
  sd.control <- apply(X = object@data[genes.use, cells.control, drop = F], 
                      MARGIN = 1, FUN = function(x) sd(x = expm1(x = x)))
  
  pooled.sd <- sqrt((((cells.control.n-1)*(sd.control^2))+((cells.treated.n-1)*(sd.treated^2)))/(cells.control.n + cells.treated.n - 2))
  
  report <- as.data.frame(matrix(nrow = length(rownames(de.test.result)), ncol = 9))
  rownames(report) <- rownames(de.test.result)
  report[,1:2] <- de.test.result[,3:4]
  report[,3] <- p.adjust(de.test.result[,1], method = "fdr")
  report[,4] <- abs((mean.treated - mean.control)/(pooled.sd))
  report[,5] <- mean.control
  report[,6] <- mean.treated
  report[,7] <- (report[,5]+0.01)/(report[,6]+0.01)
  report[,8] <- (report[,6]+0.01)/(report[,5]+0.01)
  report[,9] <- (de.test.result[,2])<0
  
  colnames(report) <- c("pct.day","pct.night","p.val.fdr","effect.size","day.avg","night.avg","day.fc","night.fc","night.up")
  temp.name <- gsub(pattern = "[+]", replacement = "", x = i)
  temp.name <- gsub(pattern = " ", replacement = "", x = temp.name)
  
  night.de.reports[[temp.name]] <- report  
  
  celltype.count = celltype.count+2
}  

night.de.reports.wilcoxon.all <- night.de.reports

#Saving and Loading Commands
while(FALSE){
save(night.de.reports.wilcoxon.all, file = "./night.de.reports.wilcoxon.all.Robj")
load("./night.de.reports.wilcoxon.all.Robj")
}

# A - Differential Expression Barplots
celltype.labels <- c("aPinealocyte","bPinealocyte","aAstrocyte","bAstrocyte","gAstrocyte","aMicroglia","bMicroglia", "VLMCs","Endothelial")

effect.size.floor <- 0.35
percent.floor <- 0.15
fold.change.floor <- 2.0
sig.cutoff <- 0.01

night.up <- c()
day.up <- c()

night.up.genes <- list()
day.up.genes <- list()

for(set in names(night.de.reports.wilcoxon.all)){
  report <- night.de.reports.wilcoxon.all[[set]]
  night.up <- c(night.up, dim(report %>% tibble::rownames_to_column('gene') %>% filter(p.val.fdr < sig.cutoff, effect.size >= effect.size.floor, pct.night >= percent.floor, 
                                                                                       night.fc >= fold.change.floor, night.up==T))[1])
  
  night.up.genes[[set]] <- (report %>% tibble::rownames_to_column('gene') %>% filter(p.val.fdr < sig.cutoff, effect.size >= effect.size.floor, pct.night >= percent.floor, 
                                                                                     night.fc >= fold.change.floor, night.up==T))$gene
  
  day.up <- c(day.up,dim(report %>% tibble::rownames_to_column('gene') %>% filter(p.val.fdr < sig.cutoff, effect.size >= effect.size.floor, pct.day >= percent.floor, 
                                                                                  day.fc >= fold.change.floor, night.up==F))[1])
  
  day.up.genes[[set]] <- (report %>% tibble::rownames_to_column('gene') %>% filter(p.val.fdr < sig.cutoff, effect.size >= effect.size.floor, pct.day >= percent.floor, 
                                                                                   day.fc >= fold.change.floor, night.up==F))$gene
}

#Aanat correction
{
  for(set in names(night.up.genes)){
    print("Aanat" %in% night.up.genes[[set]])
  }
  for(set in names(night.up.genes)[c(3:6,8)]){
    night.up.genes[[set]] <- setdiff(night.up.genes[[set]],"Aanat")
  }
  for(i in c(3:6,8)){
    night.up[i] <- night.up[i]-1
  }
}

de.bar <- data.frame(celltype=rep(celltype.labels,2), number=(as.numeric(c(night.up, day.up))))
de.bar$celltype <- factor(de.bar$celltype, levels = rev(celltype.labels))
de.bar$number[is.na(de.bar$number)] <- 0
de.bar$condition <- factor(c(rep("Night",9),rep("Day",9)), levels=c("Night","Day"))
bar.labels <- data.frame(label = c(night.up,day.up), x=(de.bar$number)+15, y=c(9:1,9:1), condition = c(rep("Night",9),rep("Day",9)))

p1 <- ggplot(data=de.bar, aes(y=number, x=celltype)) + geom_bar(stat="identity") + coord_flip() + ylab("Number of Differentially Expressed Genes") + scale_fill_manual(values=c("#424141","#A2A2A2"))
p2 <- p1 + facet_grid(. ~ condition) + theme(legend.position = "bottom", strip.background = element_blank()) + ylab("Number of Differentially Expressed Genes")
de_barplot <- p2 + geom_text(data = bar.labels, aes(x = y, y = x, label = label))

# B - DE Heatmap
genes.use <- unique(union(unlist(night.up.genes),unlist(day.up.genes)))
sig.matrix <- matrix(ncol=9, nrow=length(genes.use))
rownames(sig.matrix) <- genes.use
colnames(sig.matrix) <- celltype.labels

for(i in 1:length(names(night.de.reports.wilcoxon.all))){
  sig.genes <- night.de.reports.wilcoxon.all[[names(night.de.reports.wilcoxon.all)[i]]] %>% tibble::rownames_to_column("gene") %>% filter(night.up==T, night.fc>=2.0, effect.size>=0.35, pct.night>=0.15,p.val.fdr<0.01) %>% tibble::column_to_rownames("gene")
  sig.matrix[,i] <- rownames(sig.matrix) %in% rownames(sig.genes)
}

sig.matrix.night <- sig.matrix
sig.matrix.night[isTRUE(sig.matrix.night)] <- 1

for(i in 1:length(names(night.de.reports.wilcoxon.all))){
  sig.genes <- night.de.reports.wilcoxon.all[[names(night.de.reports.wilcoxon.all)[i]]] %>% tibble::rownames_to_column("gene") %>% filter(night.up==F, day.fc>=2.0, effect.size>=0.35, pct.day>=0.15,p.val.fdr<0.01) %>% tibble::column_to_rownames("gene")
  sig.matrix[,i] <- rownames(sig.matrix) %in% rownames(sig.genes)
}

sig.matrix.day <- sig.matrix
sig.matrix.day[isTRUE(sig.matrix.day)] <- 1
sig.matrix.day <- sig.matrix.day*-1
sig.matrix <- sig.matrix.night + sig.matrix.day
sig.matrix[is.na(sig.matrix)] <- 0

#Aanat correction
sig.matrix["Aanat",] <- c(1,1,0,0,0,0,0,0,0)
clusters <- hclust(dist(sig.matrix))
clusters <- cutree(clusters, 75)
clusters <- sort(clusters)
sig.matrix.sub <- sig.matrix[names(clusters),]
sig.matrix.m <- reshape2::melt(sig.matrix.sub)
colnames(sig.matrix.m) <- c("gene","celltype","value")

#flip cellnames
sig.matrix.m$celltype <- factor(sig.matrix.m$celltype,levels(sig.matrix.m$celltype)[9:1])

de_heatmap <- ggplot(sig.matrix.m, aes(x = gene, y = celltype)) + geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "#3066BE", mid= "light grey", high = "#D81139") +
  ylab("")
#without genes names
de_heatmap + theme(axis.text.y = element_text(angle = 45, hjust = 1), axis.text.x = element_blank())

# C - Pinealocyte Venn Diagram (These figures were manually recreated in Illustrator)
grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(x = list(night.up.genes$bPinealocyte, night.up.genes$aPinealocyte), filename = NULL, category.names = c("bPinealocyte","aPinealocyte")))
grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(x = list(day.up.genes$bPinealocyte, day.up.genes$aPinealocyte), filename = NULL, category.names = c("bPinealocyte","aPinealocyte")))

# D - Subset DE Genes Heatmap
sig.matrix.sub <- sig.matrix.sub[unique(union(unlist(night.up.genes[3:9]),unlist(day.up.genes[3:9]))),]
genes.sub <- rownames(sig.matrix.sub[rowSums(abs(sig.matrix.sub))>1,])
genes.sub <- setdiff(genes.sub, "Aanat") #Aanat correction
sig.matrix.sub <- sig.matrix.sub[genes.sub,]

clusters <- hclust(dist(sig.matrix.sub))
clusters <- cutree(clusters, 10)
clusters <- sort(clusters)

sig.matrix.sub <- sig.matrix.sub[names(clusters),]
sig.matrix.m <- reshape2::melt(sig.matrix.sub)
colnames(sig.matrix.m) <- c("gene","celltype","value")

#flip cellnames
sig.matrix.m$celltype = factor(sig.matrix.m$celltype,levels(sig.matrix.m$celltype)[9:1])

de_heatmap_sub <- ggplot(sig.matrix.m, aes(x = gene, y = celltype)) + geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "#3066BE", mid= "light grey", high = "#D81139") +
  ylab("")

de_heatmap_sub + theme(axis.text.y = element_text(angle = 45, hjust = 1), axis.text.x = element_text(angle = 45, hjust = 1, size = 6))


####
#Ouput
pdf(file = "./de_barplot.pdf", width = 3.54, height = 2.36)
de_barplot
dev.off()

#svglite(file="./de_heatmap.svg", width = 4.33, height = 2.36, pointsize = 8)
pdf(file = "./de_heatmap.pdf", width =4.33, height =2.36 , pointsize = 8)
de_heatmap + theme(axis.text.y = element_text(angle = 45, hjust = 1), axis.text.x = element_blank())
dev.off()

#svglite(file="./de_heatmap_sub.svg", width = 4.33, height = 2.36, pointsize = 8)
pdf(file = "./de_heatmap_sub.pdf", width =4.33, height =2.36 , pointsize = 8)
de_heatmap_sub
dev.off()