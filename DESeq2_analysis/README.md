# DESeq2 analysis
This document describe the use of DESeq2_analysis.R program. DESeq2 analysis is described in this paper : Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8.

## How to prepare the data ?
For each comparison, you need to prepare two files, one with the data and another with the metadata that will be stored in different folders.
### Data files
The data files need to have the following structure:  
<img width="929" alt="data" src="https://github.com/lucasDNS9/Ribes_lab/assets/127426611/070c3667-ec54-481c-afea-fdcca9194589">

The first column is containing the gene symbols, the next ones are containing the read counts for the different samples from the two conditions you want to compare. It's important to store the file as a a .txt file (

### Metadata files
The metadata files need to have the following structure:  
<img width="242" alt="metadata" src="https://github.com/lucasDNS9/Ribes_lab/assets/127426611/21a6c254-3bfe-491d-9199-364966c9d993">  

The first column contains the labels of the different samples and the second column their associated group (the two groups you want to compare). The *reference group* (the control for exemple) need to be the first group provided in the metadata file.

## How to run the code ?
The data files and metadata files need to be in the same order in their respective folder such as:  
**data_folder**  
> data_1.txt   
> data_2.txt  
> data_3.txt
  
**metadata_folder**  
> metadata_1.txt   
> metadata_2.txt  
> metadata_3.txt

**Then you can run this command line from the terminal:**
```
Rscript DESeq2_analysis.R /path/to/data_folder /path/to/metadata_folder
```
## How to interpret the results ?
You will get a .csv file as a results. The name of the file indicate the order of the comparaison: group_2_vs_group_ref.csv means that the group_2 is compared to the group_ref that is consired as the reference group (or control). The fold change is calculated as follow:  
  
$$\text{fold_change}=\log^2{(\text{group_2}\over \text{group_ref}}$$

