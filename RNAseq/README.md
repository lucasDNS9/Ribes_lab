# RNAseq analysis
This document explain how to use the python code for RNAseq analysis. It's indicated when it's necesseray to adapte the code to the specificity of your dataset.
## dataset.py
#### Section 1: import
- define working directory (folder containing code files)
- import pandas package, and filtering functions (from filtering.py file)
- import the dataset containing RNAseq data. This dataset contains for each genes the fpkm values of each conditions, and the results of the differentitial analysis (DESeq2 analysis): fold-change and adjusted p-value).
#### Section 2: Dataset organisation
#### Section 3: Row Filtering
