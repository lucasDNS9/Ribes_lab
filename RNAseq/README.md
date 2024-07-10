# RNAseq analysis
This document explain how to use the python code for RNAseq analysis. It's indicated when it's necesseray to adapte the code to the specificity of your dataset.
## dataset.py
#### Section 1: import
- define working directory (folder containing code files)
- import pandas package, and filtering functions (from filtering.py file)
- import the dataset containing RNAseq data. This dataset contains for each genes the fpkm values of each conditions, and the results of the differentitial analysis (DESeq2 analysis): fold-change and adjusted p-value).
#### Section 2: Dataset organisation
This section is used to create subset of the main dataset and keep the fpkm and comparisons in specific conditions (coolumns filtering), for exemple only with day 7 data, or only the data in BMP+ conditions. This part need to be adjusted for the specificities of your dataset.
#### Section 3: Row Filtering
This section called the functions define in filtering.py file to filter on fpkm values and extract the significantly differnetially expressed genes (DEGs).
- fpkm filtering : keep the genes (rows) if at least one of the fpkm values across conditions is above a specified threshold (here 1)
- p-value and fold-change filtering : keep the genes (rows) if the adjusted p-value is below a specified threshold (here 0.01) and if the fold-change is above a specified threshold (here 1.5). This step extract the differentially express genes (DEGs) in selected conditions.
The intermediate dataset are saved, the ones with with only fpkm filtering are used for PCA plots (cf PCA.py) or volcano plots (cf volcano.py), the DEGs datasets are used to draw heatmaps.

## filtering.py
This code define the filtering functions used to filter the rows of RNAseq dataset.
#### filter_fpkm(data, fpkm_threshold, subset='average_fpkm_')
Function to filter the rows on a minimal fpkm value. It keeps the rows if at least one of the fpkm values from specified columns is greater than the specified fpkm threshold. Arguments:
- data: dataset containing the genes, fpkm values, and DESeq results (p-values and fold-changes)
- fpkm_threshold: minimal fpkm value (usually 1)
- subset: string of character contained in the name of the fpkm columns to check. Here, I selected 'average_fpkm_' as default (average per conditions). It can be change to 'fpkm_' to select fpkm columns of every samples.
#### filter_pvalue(data, pvalue_threshold, subset='padj')
Function to filter the rows on a maximal p-value. It keeps the rows if at least one of the p-values from specified columns is below the p-value threshold. This function is not use in the main code. Arguments:
- data: dataset containing the genes, fpkm values, and DESeq results (p-values and fold-changes).
- pvalue_threshold: maximale p-value (for example 0.01)
- subset: string of character contained in the name of the p-values columns to check. Here, I selected 'padj' as default (all columns containins the adjusted p-values).
#### filter_fc(data, fc_threshold, subset='log2fc')
Function to filter the rows on a minimal fold-change. It keeps the rows if at least one of the fold-change from specified columns is greater than the fold-change threshold or lower than ...

#### filter_pvalue_fc function

## PCA.py




