# PinealGland_SingleCell
Codebook for Mays et al., 2018 (PLOS ONE). 

# Study 
Single-cell RNA sequencing of the mammalian pineal gland identifies two pinealocyte subtypes and cell type-specific daily patterns of gene expression

[https://doi.org/10.1371/journal.pone.0205883](https://doi.org/10.1371/journal.pone.0205883)

Mays JC, Kelly MC, Coon SL, Holtzclaw L, Rath MF, et al. (2018) Single-cell RNA sequencing of the mammalian pineal gland identifies two pinealocyte subtypes and cell type-specific daily patterns of gene expression. PLOS ONE 13(10): e0205883.

# Running the Analysis Code
Single-Cell RNA-Seq analysis was performed using [Seurat](github.com/satijalab/seurat) v2.2.0. 
Code should be run in the follow order:
## Setup
`SetupDependencies.R` loads libraries, `MetadataInit.R` and `GenerateMarkerReport.R` load helper functions used in subsequent scripts.

```markdown
./Scripts/SetupDependencies.R
./Scripts/MetadataInit.R
./Scripts/GenerateMarkerReport.R
```
## Single Samples
These R Notebooks contain the analysis code for each individual sample.
```
./Day1.Rmd
./Day2.Rmd
./Night1.Rmd
./Night2.Rmd
```
## Combined Samples
These R Notebooks contain the code that combines the treatment samples and applies a refined analysis. 
```
./Day1n2_NoDoublets.Rmd
./Night1n2_NoDoublets.Rmd
```
## Reference Genome Correction
The `3PrimeExtensionCorrection.R` script applies corrected count values to the Seurat objects. Corrected count values were generated by re-aligning the sequencing data using a manually updated genomic reference. 
```markdown
./Scripts/3PrimeExtensionCorrection.R
```

## Figures

These scripts were used to generate the figures for the paper.

```
./Scripts/Figure_01.R
./Scripts/Figure_03.R
./Scripts/Figure_04.R
```

