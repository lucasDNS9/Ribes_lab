# ATACseq analysis
This document explain how to use the python code for RNAseq analysis. It's indicated when it's necesseray to adapte the code to the specificity of your dataset. The datasets used in the analysis are available on the laboratory storage.
## dataset.py
#### Section 1: import
- import pandas package, and filtering functions (from filtering.py file)
- define working directory (folder containing code files)
- import the dataset containing ATACseq data. This dataset contains for each genes the rpkm values of each conditions, and the results of the differentitial analysis (DESeq2 analysis): fold-change and adjusted p-value).
#### Section 2: Dataset organisation
This section is used to create subset of the main dataset and keep the rpkm and comparisons in specific conditions (columns filtering), for exemple only with day 7 data, only the data in BMP+ conditions, or only the time-comparisons (wild-type iPSC, day 4 and day 7). This part need to be adjusted for the specificities of your dataset.
#### Section 3: Row Filtering
This section called the functions define in filtering.py file to filter on rpkm values and extract the significantly open chromatine regions (DORs).
- rpkm filtering : keep the genes (rows) if at least one of the rpkm values across conditions is above a specified threshold (here 1)
- p-value and fold-change filtering : keep the genes (rows) if the adjusted p-value is below a specified threshold (here 0.01) and if the fold-change is above a specified threshold (here 1.5). This step extract the differentially open regions (DORs) in selected conditions.

## filtering.py

