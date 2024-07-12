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
This section called the functions define in filtering.py file to filter on rpkm values and extract the significantly open chromatine regions (DORs). Datasets are saved at each step.
- rpkm filtering : keep the genes (rows) if at least one of the rpkm values across conditions is above a specified threshold (here 1)
- p-value and fold-change filtering : keep the genes (rows) if the adjusted p-value is below a specified threshold (here 0.01) and if the fold-change is above a specified threshold (here 1.5). This step extract the differentially open regions (DORs) in selected conditions.

## filtering.py
This code define the filtering functions used to filter the rows of ATACseq dataset.
#### filter_rpkm(data, rpkm_threshold, subset='average_RPKM_')
Function to filter the rows on a minimal rpkm value. It keeps the rows if at least one of the rpkm values from specified columns is greater than the specified rpkm threshold. \
Arguments:
- data: dataset containing the genes, rpkm values, and DESeq results (p-values and fold-changes)
- rpkm_threshold: minimal rpkm value (usually 1)
- subset: string of character contained in the name of the rpkm columns to check. Here, I selected 'average_RPKM_' as default (average per conditions). It can be change to 'RPKM_' to select rpkm columns of every samples.
#### filter_pvalue(data, pvalue_threshold, subset='padj')
Function to filter the rows on a maximal p-value. It keeps the rows if at least one of the p-values from specified columns is below the p-value threshold. This function is not use in the main code. \
Arguments:
- data: dataset containing the genes, fpkm values, and DESeq results (p-values and fold-changes)
- pvalue_threshold: maximale p-value (for example 0.01)
- subset: string of character contained in the name of the p-values columns to check. Here, I selected 'padj' as default (all columns containins the adjusted p-values).
#### filter_fc(data, fc_threshold, subset='log2fc')
Function to filter the rows on a minimal fold-change. It keeps the rows if at least one of the fold-change from specified columns is greater than the fold-change threshold or lower than 1/{fold-change threshold}. This function is not use in the main code. \
Arguments:
- data: dataset containing the genes, rpkm values, and DESeq results (p-values and fold-changes)
- fc_threshold: fold-change threshold (usually 1.5)
- subset: string of character contained in the name of the fold-change columns to check. Here, I selected 'log2fc' as default (all columns containins the log2 of fold-changes)\
Special Warning: the DESeq analysis gives fold-changes as $log2(fold change)$, this specificity is taken into account in the function. If it's not the case in your dataset, you need to adjuste this function.
#### filter_pvalue_fc(data, pvalue_threshold, fc_threshold, subset_p='padj', subset_fc='log2fc')
This function filter the rows if both the p-value and fold-change of one of the comparisons meet the selection criteria. This function is better to extract the differentially expressed genes compared to the successive use of the two previous functions.\
Arguments:
- data: dataset containing the genes, rpkm values, and DESeq results (p-values and fold-changes)
- pvalue_threshold: maximale p-value (for example 0.01)
- fc_threshold: fold-change threshold (usually 1.5), keep the row if at least one of the fold-change from specified columns is greater than the fold-change threshold or lower than 1/{fold-change threshold}. It's important to notice that the DESeq analysis gives fold-changes as $log2(fold change)$, this specificity is taken into account in the function. If it's not the case in your dataset, you need to adjuste this function.
- subset_p: string of character contained in the name of the p-values columns to check. Here, I selected 'padj' as default (all columns containins the adjusted p-values)
- subset_fc: string of character contained in the name of the fold-change columns to check. Here, I selected 'log2fc' as default (all columns containins the log2 of fold-changes)\
Specificity: The function will check both p-values and fold-change at the same time for each rows. If there is several comparisons to check, the p-value columns and fold-change columns need to be in the same order : comparison_1 > comparison_2 > comparison_3. Cf. table organization presented above.

## regions_to_genes.py
#### Section 1: Import
- Import packages
- Define working directory
- Import the table with the associations between a chroomatin regions and the two closest genes. Organization of this table with a column 'chrom:start-end' as the row indexes.
- Import the datasets containings the Differentially Open Regions (DORs) extract with the filtering functions.
#### Section 2: Define the function
This function associate the chromatin regions (from the DORs list) and the two closest genes using the association table. It generate a new dataset with the region locations, the rpkm values and DESeq results (p-value and fold-change), and the closests genes. Function: **region_to_gene(data_association, data_DORs)**. \
 - data_association: dataset associating chromatin region and two closest genes (2closest_atac.bed)
 - data_DORs: datasets containing the location of the DORs
#### Section 3: Application of the function
Apply the function the the DORs from the differents comparisons and save the results as new tables.

## timecomp.py
This code is used to represente the chromatin opening modification during development, comparing the ATACseq signal (rpkm) between different timepoints.
#### Section 1: Import
- Import required packages
- Define working directory
- Import datasets containing the Differentially Open Regions (DORs) between timepoints (Wild-types between iPSCs at day 0 and organoids at day 4 and day 7). These datasets are generated with the code dataset.py and filtering.py (cf. above)
- Define the list of columns (their names) containing the p-values of the different time-comparisons and fold-change of the different time-comparisons (they have to be in the same order: comparison_1 > comparison_2 > comparison_3). This has to be adapted to the specificities of your dataset.
#### Section 2: Classification function
The function classifies each DORs into one 'profil of opening' according to the results of DESeq analysis between the different timepoint comparisons. For each comparisons, the function check wheter theiir is a significant difference (based on a p-value threshold and a fold-change threshold) and associate a letter: 'o'=more opened compared to the reference group, 'c'=more closed compared to the reference group, and 'n'=no significant differences. The reference group is define during DESeq analysis (cf.corresponding folder on GitHub). Function: **classification(data, p_comparisons, fc_comparisons, p=0.01, fc=1.5)** \
Arguments:
- data: dataset containing the chromatin regions ('chr:start-end') and the DESeq results (p-values and fold-change) for each time comparisons. The function will create a new columns to this dataset containing the 'opening profils'
- p_comparisons: list of the columns containing the p-value of the different time-comparisons
- fc_comparisons: list of the columns containing the fold-change of the different time-comparisons
- p: p-value threshold of significance (significant if p-value < threshold), the default value is 0.01.
- fc: fold-change threshold (significant if fold-change > threshold or < -threshold), the default value is 1.5. The DESeq results give the $log2(fold change)$, it's taken into account in the function.
#### Section 3: Analysis and processing
- Print the number of regions for each profils
- Remove the aberrant profils from the classification, for exemple: '1o2n3c','1o2o3n', and '1c2n3o'
- Process the data for graphical representation: set the 'chr:start-end' column as the row indexes, only keep the columns with the average rpkm and the profils
#### Section 4: Graphs
Define the graphical function to draw lineplots for each profil, representing the median±sd rpkm. The function save automatically the graphs in the working directory (here 18 graphs). Function: **plotting_profils(DORs_plot)**, the argument DORs_plot is the dataset generated in **section 3** containing the chromatin region as the row indexes, the average rpkm for the different timepoint (here day 0, day 4 and day 7), and a column with the opening profils.

## highlight_timeplot.py
This code is used to highlight selected chromatin regions in the profil plots representing the chromatin opening changes during development (comparisons of ATACseq signal between different timepoints). Here, i want to highlight the chromatin regions nearby genes identified in the RNAseq analysis, but this code can be used to highlight any chromatin regions from a given list in the profil plots.
#### Section 1: Import
- Import the required packages
- Define the working directory
- Import the dataset containing the chromatin regions differentially open (DORs) during development (time-comparisons) and set the chromatin region location as the row indexes.
- Import the DORs between homozygous and wild-type in BMP conditions (Here, I want to highlight in these regions the one nearby Differentially Expressed Genes (DEGs) identified in RNAseq analysis.
- Import the DEGs from RNAseq analysis
- Subset the chromatin regions nearby the genes contains in the list of DEGs
This section can be modify to extract among all the genes nearby chromatin region dynamic during development genes that are in a specified list.

#### Section 2: Plotting function
This code is used to represente the chromatin opening modification during development, comparing the ATACseq signal (rpkm) between different timepoints, and highlight the chromatin regions nearby specified genes (here, the genes that were identified in the RNAseq analysis). Function: **plotting_profils_list(DORs_plot, list_regions)**. \
Arguments:
- DORs_plot: dataset containing the chromatin region location (row indexes), the average rpkm for the different timepoints, and the opening profils. 
- list_regions: list of chromatin regions to highlight
This function can be adjusted according to the specificities of your data (different timepoints ...)




