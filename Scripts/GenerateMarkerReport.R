GenerateMarkerReport <- function(object = NULL, ident.1 = NULL, ident.2 = NULL, marker.report.raw = NULL){
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
  marker.report$pct.2 <- marker.report.raw[,4]
  marker.report$p.val.fdr <- p.adjust(marker.report.raw[,1], method = "fdr")
  marker.report$effect.size <- abs((mean.ident.2 - mean.ident.1)/(pooled.sd))
  marker.report$mean.1 <- mean.ident.1
  marker.report$mean.2 <- mean.ident.2
  marker.report$log.fc <- marker.report.raw[,2]
  marker.report$ident.1.fc <- (marker.report$mean.1+0.01)/(marker.report$mean.2+0.01)
  marker.report$ident.2.fc <- (marker.report$mean.2+0.01)/(marker.report$mean.1+0.01)
  marker.report$ident.2.up <- (marker.report.raw[,2])<0
  marker.report$gene <- rownames(marker.report.raw)

  return(marker.report)
}