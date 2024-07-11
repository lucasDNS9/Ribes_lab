# RNAseq analysis
This document explain how to use the python code for RNAseq analysis. It's indicated when it's necesseray to adapte the code to the specificity of your dataset. The datasets used in the analysis are available on the laboratory storage.
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
Function to filter the rows on a minimal fpkm value. It keeps the rows if at least one of the fpkm values from specified columns is greater than the specified fpkm threshold. \
Arguments:
- data: dataset containing the genes, fpkm values, and DESeq results (p-values and fold-changes)
- fpkm_threshold: minimal fpkm value (usually 1)
- subset: string of character contained in the name of the fpkm columns to check. Here, I selected 'average_fpkm_' as default (average per conditions). It can be change to 'fpkm_' to select fpkm columns of every samples.
#### filter_pvalue(data, pvalue_threshold, subset='padj')
Function to filter the rows on a maximal p-value. It keeps the rows if at least one of the p-values from specified columns is below the p-value threshold. This function is not use in the main code. \
Arguments:
- data: dataset containing the genes, fpkm values, and DESeq results (p-values and fold-changes)
- pvalue_threshold: maximale p-value (for example 0.01)
- subset: string of character contained in the name of the p-values columns to check. Here, I selected 'padj' as default (all columns containins the adjusted p-values).
#### filter_fc(data, fc_threshold, subset='log2fc')
Function to filter the rows on a minimal fold-change. It keeps the rows if at least one of the fold-change from specified columns is greater than the fold-change threshold or lower than 1/{fold-change threshold}. This function is not use in the main code. \
Arguments:
- data: dataset containing the genes, fpkm values, and DESeq results (p-values and fold-changes)
- fc_threshold: fold-change threshold (usually 1.5)
- subset: string of character contained in the name of the fold-change columns to check. Here, I selected 'log2fc' as default (all columns containins the log2 of fold-changes)\
Special Warning: the DESeq analysis gives fold-changes as $log2(fold change)$, this specificity is taken into account in the function. If it's not the case in your dataset, you need to adjuste this function.
#### filter_pvalue_fc(data, pvalue_threshold, fc_threshold, subset_p='padj', subset_fc='log2fc')
This function filter the rows if both the p-value and fold-change of one of the comparisons meet the selection criteria. This function is better to extract the differentially expressed genes compared to the successive use of the two previous functions.\
Arguments:
- data: dataset containing the genes, fpkm values, and DESeq results (p-values and fold-changes)
- pvalue_threshold: maximale p-value (for example 0.01)
- fc_threshold: fold-change threshold (usually 1.5), keep the row if at least one of the fold-change from specified columns is greater than the fold-change threshold or lower than 1/{fold-change threshold}. It's important to notice that the DESeq analysis gives fold-changes as $log2(fold change)$, this specificity is taken into account in the function. If it's not the case in your dataset, you need to adjuste this function.
- subset_p: string of character contained in the name of the p-values columns to check. Here, I selected 'padj' as default (all columns containins the adjusted p-values)
- subset_fc: string of character contained in the name of the fold-change columns to check. Here, I selected 'log2fc' as default (all columns containins the log2 of fold-changes)\
Specificity: The function will check both p-values and fold-change at the same time for each rows. If there is several comparisons to check, the p-value columns and fold-change columns need to be in the same order : comparison_1 > comparison_2 > comparison_3. Cf. table organization presented above.

## PCA.py
This code was used to do a Principal Component Analysis (PCA) from RNAseq data and draw the following graph:
[PCA_d7.pdf](https://github.com/user-attachments/files/16163061/PCA_d7.pdf)
#### Section 1: Data preprocessing
- Import required packages
- Define working directory (folder containing the datasets to import)
- Import the datasets containing the data filtered on fpkm values (previously saved from dataset.py), namely the whole dataset and only the data at day 7 (represented in the graph above). The datasets contains a column with gene names and columns with fpkm for each conditions as follow:
  
<img width="1081" alt="Capture d’écran 2024-07-10 à 15 59 48" src="https://github.com/lucasDNS9/Ribes_lab/assets/127426611/6f453862-2d82-4135-8668-499e255ad54f"> 

- Import a dataset containing the association between samples and condition groups, it's important to keep the structure as follow:
  
<img width="299" alt="Capture d’écran 2024-07-10 à 16 06 15" src="https://github.com/lucasDNS9/Ribes_lab/assets/127426611/6e0d6739-c8f9-4972-a087-6bd14290e966">

- Define the gene names as the row indexes.
- Logarithmic transformationof the fpkm values: improve the graphical representation
#### Section 2: Graphical functions
- Import packages (pyplot from matplotlib and decomposition from scikit-learn)
- Define **scree_plot** function. This function draw a scree plot (contribution of the different composant from the PCA), that look like that: [scree_plot_d7.pdf](https://github.com/user-attachments/files/16163671/scree_plot_d7.pdf). function: **scree_plot(data, subset, title=None)**:
  - data: dataset containing fpkm values (filtered on minimal fpkm)
  - subset: string of character contained in the name of the fpkm columns (usually 'fpkm_')
  - title: if a title is specified such as title='title.pdf', the graph will be saved with the specified title (title=None by default).
- Define PCA function, **plot_PCA(data, subset, sample_group, n_components=2, title=None)**:
  - data: dataset containing fpkm values (filtered on minimal fpkm)
  - subset: string of character contained in the name of the fpkm columns (usually 'fpkm_')
  - sample_groupe: dataset containing the association between samples and condition group (previously imported as sample_groupe)
  - n_components: number of principal components (PC) to include in the PCA, the default value is 2, the graph will always display 2 components but you can change to higher number to have in addition to the graph a table with the contribution more PC for each samples.
  - title: if a title is specified such as title='title.pdf', the graph will be saved with the specified title (title=None by default)
In the PCA function is defined two dictionnaries that are used to improve graphic epresentation: color_mapping that map condition groups with a color (here the genotype), and shape_mapping that map condition groups to a shape (here the presence of abscence of BMP). Those two dictionnaries and the way they map the different conditions to color and shape need to be adjusted to the specificities of your dataset.
#### Section 3: Graphs
Application of the functions scree_plot() and plot_PCA().

## volcano.py
This code was used to draw volcano plots from RNAseq data.
#### Section 1: Data preprocessing
- import required packages
- define working directory
- Import the datasets containing the data filtered on fpkm values (previously saved from dataset.py), with gene names, fpkm for each condition, p-values and fold-changes for the different comparisons. Here two datasets : the gene regulated by PAX3 (PAX3 DEGs) in the two conditions with and without BMP.
- Define the gene names as the row indexes.
#### Section 2: Functions
- Define **select_comp()** function. This function extract the data required to draw the volcano plot for only one comparison, namely the columns containing the adjusted p-value ('padj') and the fold-change ('log2fc) of the specified comparisons, this can be change if necessary. This function is use in the volcano plot function. Function: **select_comp(data, condition_1, condition_2)**:
  - data: dataset containing fpkm values (filtered on minimal fpkm)
  - condition_1: the first condition we want to compare
  - condition_2: the second condition we want to compare (the order between condition 1 and 2 does not matter, the order of the comparison during DESeq analysis.
- Define **gene_extraction()** function. This function extract the Differentially Expressed Genes (DEGs) from a single comparison. It can be used separetely from the other functions but it's used in the volcano plot function. Function: **gene_extraction(data, condition_1, condition_2, fc_threshold, pvalue_threshold, save=False** :
  - data: dataset containing the p-values and fold-changes for one comparisons (the gene names are set as the row indexes), this table is generate with the previous **select_comp()** function, it's taken into account so the dataset can contains other information before using this function
  - conditions_1: the first condition we want to compare
  - condition_2: the second condition we want to compare (the order between condition 1 and 2 does not matter)
  - fc_threshold: fold-change threshold (usually 1.5), the DESeq analysis gives fold-changes as $log2(fold change)$, this specificity is taken into account in the function. If it's not the case in your dataset, you need to adjuste this function.
  - pvalue_threshold: maximale p-value (for example 0.01)
  - save: if save=True, a .csv file is saved with a list of DEGs and their associated p-values and fold-changes from the specified comparison.
- Define **volcano_plot()** function. This function is used to represent the p-value and fold-change of each genes and highlight the one that are differentially expressed between two conditions, for example: [volcano.pdf](https://github.com/user-attachments/files/16179792/hom_BMP_D7_vs_ctl_BMP_D7.pdf). Function: **volcano_plot(data, condition_1, condition_2, fc_threshold, pvalue_threshold=0.01, save=False, pdf=False)**:
  - data: dataset containing the p-values and fold-changes for one comparisons (the gene names are set as the row indexes), this table is generate with the previous **select_comp()** function, it's taken into account so the dataset can contains other information before using this function
  - condition_1: the first condition we want to compare
  - condition_2: the second condition we want to compare (the order between condition 1 and 2 does not matter)
  - fc_threshold: fold-change threshold (usually 1.5), the DESeq analysis gives fold-changes as $log2(fold change)$, this specificity is taken into account in the function. If it's not the case in your dataset, you need to adjuste this function.
  - pvalue-threshold: maximale p-value (for example 0.01)
  - save: if save=True, it saves a pdf with the gene differentially expressed (DEGs, highlighted in the graph).
  - pdf: if pdf=True, it saves a pdf with the graph.
#### Section 3: Generate the Volcano plots
Exemples of how the function has beed used to generate the volcano plots


