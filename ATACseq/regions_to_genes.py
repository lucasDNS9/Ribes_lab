#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:59:40 2024

@author: lucasdenis
"""
'''
Section 1: Import
'''
#import packages
import os
import pandas as pd

#Define working directory
os.chdir('/Users/lucasdenis/Desktop/ATACseq')

#import table with associations regions-genes
close_genes = pd.read_csv('./data/2closest_atac.bed', sep='\t', header=None)
#set up column index
close_genes.columns = close_genes.iloc[0]
close_genes = close_genes.drop(close_genes.index[0])
#fuse the columns 'chrom', 'start', 'end' and set as row index
close_genes['chrom:start-end'] = close_genes['chrom'] + ':' + close_genes['start'].astype(str) + '-' + close_genes['end'].astype(str)
close_genes.set_index('chrom:start-end', inplace=True)

#import DORs datasets
ATACseq_PAX3_BMP_DORs = pd.read_csv('./results/lists_DORs/ATACseq_PAX3_BMP_DORs.csv', sep=';', decimal='.', index_col=0)
ATACseq_PAX3_NO_DORs = pd.read_csv('./results/lists_DORs/ATACseq_PAX3_NO_DORs.csv', sep=';', decimal='.', index_col=0)
ATACseq_timecomp_DORs = pd.read_csv('./results/lists_DORs/ATACseq_timecomp_DORs.csv', sep=';', decimal='.', index_col=0)
ATACseq_hom_vs_WT_BMP_DORs =  pd.read_csv('./results/lists_DORs/ATAseq_hom_vs_WT_BMP.csv', sep=';', decimal='.', index_col=0)

###############################################################################
'''
Section 2: Define the function
'''
def region_to_gene(data_association, data_DORs):
    
    #cross the DORs with the association dataframe
    data = data_association[(data_association.index).isin(data_DORs.index)]
    
    #repeat twice the DORs dataframe (each region associated with two genes)
    data_DORs_2 = data_DORs.loc[data_DORs.index.repeat(2)].reset_index(drop=True)
    
    #concat dataframes and set up DORs labels
    data_DORs_genes = pd.concat([data.reset_index(drop=True),data_DORs_2], axis=1)
    data_DORs_genes['chrom:start-end'] = data_DORs_genes['chrom'] + ':' + data_DORs_genes['start'].astype(str) + '-' + data_DORs_genes['end'].astype(str)
    data_DORs_genes.set_index('chrom:start-end', inplace=True)

    return data_DORs_genes


###############################################################################
'''
Section 3: Application of the function
'''
#create the tables with the associations DORs-genes
ATACseq_PAX3_BMP_DORs_genes = region_to_gene(data_association=close_genes, data_DORs=ATACseq_PAX3_BMP_DORs)
ATACseq_PAX3_NO_DORs_genes = region_to_gene(data_association=close_genes, data_DORs=ATACseq_PAX3_NO_DORs)
ATACseq_timecomp_DORs_genes = region_to_gene(data_association=close_genes, data_DORs=ATACseq_timecomp_DORs)
ATACseq_hom_vs_WT_DORs_genes = region_to_gene(data_association=close_genes, data_DORs=ATACseq_hom_vs_WT_BMP_DORs)

#Save the results
ATACseq_PAX3_BMP_DORs_genes.to_csv('ATACseq_PAX3_BMP_DORs_genes.csv')
ATACseq_PAX3_NO_DORs_genes.to_csv('ATACseq_PAX3_NO_DORs_genes.csv')
ATACseq_timecomp_DORs_genes.to_csv('ATACseq_timecomp_DORs_genes.csv')
ATACseq_hom_vs_WT_DORs_genes.to_csv('ATACseq_hom_vs_WT_BMP_DORs_genes.csv')




