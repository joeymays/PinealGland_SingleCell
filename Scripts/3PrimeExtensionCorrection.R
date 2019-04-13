edited.genes <- c("Chrnb2","Clic4","Crx","Gabrg2","Hcn2","Kcnab2","Kcnq4","Map2k4","Map2k7","Map4k5","Mapk1","P2ry4","Prkaa2","Prkacb","Prkd1","Scn7a","Tab2","Adrb1","Clcn4","Gabrb3","Hcn1","Ppm1a","Ppp1r7","Prkar2a","Prkce","Rbp3","Scn2b","Scn8a")

Day1.edit.raw<-Read10X("./ExtendedReferenceRawMatrix/Day1/")
Day2.edit.raw<-Read10X("./ExtendedReferenceRawMatrix/Day2/")

Day1.edit<-CreateSeuratObject(Day1.edit.raw)
Day2.edit<-CreateSeuratObject(Day2.edit.raw)

Day1n2.edit<-MergeSeurat(Day1.edit, Day2.edit, add.cell.id1 = "D1", add.cell.id2 = "D2", do.normalize = F)

edited.genes <- edited.genes[edited.genes %in% rownames(Pineal2_Day1n2@raw.data)]

Day1n2.edit.matrix <- Day1n2.edit@raw.data[edited.genes,colnames(Pineal2_Day1n2@raw.data)]
Day1n2.edit.matrix <- as.matrix(Day1n2.edit.matrix)

for(i in edited.genes){
  Pineal2_Day1n2@raw.data[i,] <- Day1n2.edit.matrix[i,] 
  print(i)
}

Night1.edit.raw<-Read10X("./ExtendedReferenceRawMatrix/Night1/")
Night2.edit.raw<-Read10X("./ExtendedReferenceRawMatrix/Night2/")

Night1.edit<-CreateSeuratObject(Night1.edit.raw)
Night2.edit<-CreateSeuratObject(Night2.edit.raw)

Night1n2.edit<-MergeSeurat(Night1.edit, Night2.edit, add.cell.id1 = "N1", add.cell.id2 = "N2", do.normalize = F)

edited.genes <- c("Chrnb2","Clic4","Crx","Gabrg2","Hcn2","Kcnab2","Kcnq4","Map2k4","Map2k7","Map4k5","Mapk1","P2ry4","Prkaa2","Prkacb","Prkd1","Scn7a","Tab2","Adrb1","Clcn4","Gabrb3","Hcn1","Ppm1a","Ppp1r7","Prkar2a","Prkce","Rbp3","Scn2b","Scn8a")
edited.genes <- edited.genes[edited.genes %in% rownames(Pineal2_Night1n2@raw.data)]

Night1n2.edit.matrix <- Night1n2.edit@raw.data[edited.genes,colnames(Pineal2_Night1n2@raw.data)]
Night1n2.edit.matrix <- as.matrix(Night1n2.edit.matrix)

for(i in edited.genes){
  print(i)
  Pineal2_Night1n2@raw.data[i,] <- Night1n2.edit.matrix[i,] 
}

Pineal2_Day1n2 <- NormalizeData(Pineal2_Day1n2)
Pineal2_Night1n2 <- NormalizeData(Pineal2_Night1n2)

Pineal2_Combined.2 <- MergeSeurat(Pineal2_Day1n2, Pineal2_Night1n2)
celltypes <- levels(Pineal2_Day1n2@ident)
for(i in celltypes){
  Pineal2_Combined.2 <- SetIdent(Pineal2_Combined.2, cells.use = WhichCells(Pineal2_Day1n2, i), ident.use = i)
  Pineal2_Combined.2 <- SetIdent(Pineal2_Combined.2, cells.use = WhichCells(Pineal2_Night1n2, i), ident.use = i)
}
Pineal2_Combined.2@ident = factor(Pineal2_Combined.2@ident,levels(Pineal2_Combined.2@ident)[c(2,8,4,5,6,1,7,3,9)])
