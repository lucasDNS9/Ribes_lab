#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:26:18 2024

@author: lucasdenis
"""
###############################################################################
'''
set working directory, import packages, functions and dataset
'''

#define working directory (containing the filtering.py file (with filtering functions))
import os
os.chdir('/Users/lucasdenis/Desktop/RNAseq/code')

#import packages and functions
import pandas as pd
from filtering import filter_fpkm, filter_pvalue_fc

#import dataset
RNAseq_data = pd.read_csv('/Users/lucasdenis/Desktop/RNAseq/datasets/RNAseq_data.csv', sep=';', decimal='.')

###############################################################################
'''
Dataset organization and subset (column filtering)
'''

#define the dataset row indexes
RNAseq_data.set_index('gene_symbol', inplace=True)


#Define the datasets i will be working with (subsets of main dataset)

#only d7 data and associated comparisons
RNAseq_data_d7 = RNAseq_data[[col for col in RNAseq_data.columns if '_iPS' not in col and '_d4' not in col]]
#only d7 with BMP data and associated comparisons
RNAseq_BMP = RNAseq_data_d7[[col for col in RNAseq_data_d7.columns if '_B' in col]]
#only d7 without BMP data and associated comparisons
RNAseq_NO = RNAseq_data_d7[[col for col in RNAseq_data_d7.columns if '_N' in col]]

#datasets with only PAX3 effects and associated comparisons in BMP or NO conditions
RNAseq_PAX3_BMP = RNAseq_BMP[[col for col in RNAseq_BMP.columns if 'NO' not in col]]
RNAseq_PAX3_NO = RNAseq_NO[[col for col in RNAseq_NO.columns if 'BMP' not in col]]


###############################################################################
'''
Rows Filtering (1={fpkm}, 2={fpkm + pvalue + fold change})
'''
#filter on fpkm (genes with at least on fpkm > fpkm_threshold remain)

#main dataset
RNAseq_filtered = filter_fpkm(RNAseq_data, fpkm_threshold=1)

#only d7
RNAseq_d7_filtered = filter_fpkm(RNAseq_data_d7, fpkm_threshold=1)

#either BMP or NO BMP
RNAseq_PAX3_BMP_1 = filter_fpkm(RNAseq_PAX3_BMP, fpkm_threshold=1, subset='average_fpkm')
RNAseq_PAX3_NO_1 = filter_fpkm(RNAseq_PAX3_NO, fpkm_threshold=1, subset='average_fpkm')

#save the datasets
RNAseq_filtered.to_csv('RNAseq_fpkm_filt.csv')
RNAseq_d7_filtered.to_csv('RNAseq_d7_fpkm_filt.csv')
RNAseq_PAX3_BMP_1.to_csv('RNAseq_PAX3_BMP_fpkm_filt.csv')
RNAseq_PAX3_NO_1.to_csv('RNAseq_PAX3_NO_fpkm_filt.csv')

###

#filter on p-value and fold change

#Differentially Expressed Genes (DEGs) in BMP conditon
RNAseq_PAX3_BMP_2 = filter_pvalue_fc(RNAseq_PAX3_BMP_1,pvalue_threshold=0.01, fc_threshold=1.5)
#DEGs in NO BMP condition
RNAseq_PAX3_NO_2 = filter_pvalue_fc(RNAseq_PAX3_NO_1,pvalue_threshold=0.01, fc_threshold=1.5)

#save the datasets
RNAseq_PAX3_BMP_2.to_csv('RNAseq_PAX3_BMP_DEGs.csv')
RNAseq_PAX3_NO_2.to_csv('RNAseq_PAX3_NO_DEGs.csv')




