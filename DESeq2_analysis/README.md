# DESeq2 analysis
This document describe the use of DESeq2_analysis.R program. DESeq2 analysis is described in this paper : Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8.

## How to prepare the data ?
For each comparison, you need to prepare two files, one with the data and another with the metadata that will be stored in different folders.
### Data files
The data files need to have the following structure:
<img width="813" alt="Capture d’écran 2024-04-24 à 15 00 44" src="https://github.com/lucasDNS9/Ribes_lab/assets/127426611/b9f92b2d-3d04-441e-adb0-36dbc4a8bc68">
### Metadata files
The metadata files need to have the following structure:

## How to run the code ?

```
Rscript DESeq2_analysis.R /path/to/data_folder /path/to/metadata_folderD
```
