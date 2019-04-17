### 10X Genomics Initial Metadata Processing 
#Input: Seurat Object v2.3.4
#Stores BARCODES in Metadata table column 'barcode'
#Creates Cell IDs from a given prefix and numbers in increments of 1
#Stores Cell IDs in Metadata Table
### ///////////////////////////////

MetadataInit <- function(SeuratObject, ID.prefix){
  
  # Save Cell IDs to MetaData Slot
  barcodes <- colnames(GetAssayData(SeuratObject, assay = "RNA", slot = "counts"))
  SeuratObject@meta.data[,"barcode"] <- barcodes
  
  cell.ids = c()
  for(i in 1:length(colnames(GetAssayData(SeuratObject, assay = "RNA", slot = "counts")))){
    cell.id.current <- paste(ID.prefix, i, sep = '_')
    cell.ids <- c(cell.ids, cell.id.current)
  }  
  SeuratObject@meta.data[,"cell.id"] <- cell.ids
  return(SeuratObject)
}
