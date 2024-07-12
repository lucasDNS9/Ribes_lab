#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:26:18 2024

@author: lucasdenis
"""
'''
Section 1 : Import
'''
#import packages and function
import os
import pandas as pd
from filtering import filter_rpkm, filter_pvalue_fc

#Define working directory
os.chdir('/Users/lucasdenis/Desktop/ATACseq')

#import dataset
ATACseq_data = pd.read_csv('./data_ATACseq.csv', sep=';', decimal='.')

###############################################################################
'''
Section 2: Dataset organization
'''
#define the index
ATACseq_data.set_index('chr:start-end', inplace=True)


#Define the datasets I will be working with (subsets of main dataset)

#only d7 data and associated comparisons
ATACseq_data_d7 = ATACseq_data[[col for col in ATACseq_data.columns if '_iPS' not in col and '_d4' not in col]]
#only d7 with BMP data and associated comparisons
ATACseq_BMP = ATACseq_data_d7[[col for col in ATACseq_data_d7.columns if 'BMP' in col]]
#only d7 without BMP data and associated comparisons
ATACseq_NO = ATACseq_data_d7[[col for col in ATACseq_data_d7.columns if 'NO' in col]]

#datasets with only PAX3 effects and associated comparisons in BMP or NO conditions
ATACseq_PAX3_BMP = ATACseq_BMP[[col for col in ATACseq_BMP.columns if 'NO' not in col]]
ATACseq_PAX3_NO = ATACseq_NO[[col for col in ATACseq_NO.columns if 'BMP' not in col]]

#dataset with only homozygous and WT in BMP
ATACseq_hom_vs_WT_BMP = ATACseq_PAX3_BMP[[col for col in ATACseq_PAX3_BMP.columns if 'het' not in col]]
#dataset with only heterozygous and WT in BMP
ATACseq_het_vs_WT_BMP = ATACseq_PAX3_BMP[[col for col in ATACseq_PAX3_BMP.columns if 'hom' not in col]]
#dataset with only hom vs het in BMP
ATACseq_hom_vs_het_BMP = ATACseq_PAX3_BMP[[col for col in ATACseq_PAX3_BMP.columns if 'ctl' not in col]]

#dataset with only homozygous and WT in NO
ATACseq_hom_vs_WT_NO = ATACseq_PAX3_NO[[col for col in ATACseq_PAX3_NO.columns if 'het' not in col]]
#dataset with only heterozygous and WT in NO
ATACseq_het_vs_WT_NO = ATACseq_PAX3_NO[[col for col in ATACseq_PAX3_NO.columns if 'hom' not in col]]
#dataset with only hom vs het in NO
ATACseq_hom_vs_het_NO = ATACseq_PAX3_NO[[col for col in ATACseq_PAX3_NO.columns if 'ctl' not in col]]


#only time comparisons (BMP and NO)
ATACseq_timecomp_BMP = ATACseq_data[[col for col in ATACseq_data.columns if 'hom' not in col and 'het' not in col and "NO" not in col and 'BMP_d7_vs_ctl' not in col]]
ATACseq_timecomp_NO = ATACseq_data[[col for col in ATACseq_data.columns if 'hom' not in col and 'het' not in col and "BMP" not in col and 'BMP_d7_vs_ctl' not in col]]


###############################################################################
'''
Rows Filtering (1={rpkm}, 2={rpkm + pvalue + fold change})
'''
#filter on fpkm (genes with at least one average rpkm per condition > rpkm_threshold)

#main dataset
ATACseq_filtered = filter_rpkm(ATACseq_data, rpkm_threshold=1, subset='average_RPKM')

#only d7
ATACseq_d7_filtered = filter_rpkm(ATACseq_data_d7, rpkm_threshold=1, subset='average_RPKM')

#either BMP or NO BMP
ATACseq_PAX3_BMP_1 = filter_rpkm(ATACseq_PAX3_BMP, rpkm_threshold=1, subset='average_RPKM')
ATACseq_PAX3_NO_1 = filter_rpkm(ATACseq_PAX3_NO, rpkm_threshold=1, subset='average_RPKM')

#hom_vs_WT, het_vs_WT and hom_vs_het in BMP
ATACseq_hom_vs_WT_BMP_1 = filter_rpkm(ATACseq_hom_vs_WT_BMP, rpkm_threshold=1, subset='average_RPKM')
ATACseq_het_vs_WT_BMP_1 = filter_rpkm(ATACseq_het_vs_WT_BMP, rpkm_threshold=1, subset='average_RPKM')
ATACseq_hom_vs_het_BMP_1 = filter_rpkm(ATACseq_hom_vs_het_BMP, rpkm_threshold=1, subset='average_RPKM')

#hom_vs_WT, het_vs_WT and hom_vs_het in NO
ATACseq_hom_vs_WT_NO_1 = filter_rpkm(ATACseq_hom_vs_WT_NO, rpkm_threshold=1, subset='average_RPKM')
ATACseq_het_vs_WT_NO_1 = filter_rpkm(ATACseq_het_vs_WT_NO, rpkm_threshold=1, subset='average_RPKM')
ATACseq_hom_vs_het_NO_1 = filter_rpkm(ATACseq_hom_vs_het_NO, rpkm_threshold=1, subset='average_RPKM')


#time comparisons BMP OR NO
ATACseq_timecomp_BMP_1 = filter_rpkm(ATACseq_timecomp_BMP, rpkm_threshold=1, subset='average_RPKM')
ATACseq_timecomp_NO_1 = filter_rpkm(ATACseq_timecomp_NO, rpkm_threshold=1, subset='average_RPKM')

#save the datasets
ATACseq_filtered.to_csv('ATACseq_rpkm_filt.csv')
ATACseq_d7_filtered.to_csv('ATACseq_d7_rpkm_filt.csv')
ATACseq_PAX3_BMP_1.to_csv('ATACseq_PAX3_BMP_rpkm_filt.csv')
ATACseq_PAX3_NO_1.to_csv('ATACseq_PAX3_NO_rpkm_filt.csv')

###############################################################################

#filter on p-value and fold change

#Differentially Open Regions (DORs) in BMP condition
ATACseq_PAX3_BMP_2 = filter_pvalue_fc(ATACseq_PAX3_BMP_1,pvalue_threshold=0.01, fc_threshold=1.5, 
                                      subset_p='padj', subset_fc='log2fc')
#DORs in NO BMP condition
ATACseq_PAX3_NO_2 = filter_pvalue_fc(ATACseq_PAX3_NO_1,pvalue_threshold=0.01, fc_threshold=1.5,
                                     subset_p='padj', subset_fc='log2fc')

#DORs per comparison in BMP
ATACseq_hom_vs_WT_BMP_2 = filter_pvalue_fc(ATACseq_hom_vs_WT_BMP_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                           subset_p='padj', subset_fc='log2fc')
ATACseq_het_vs_WT_BMP_2 = filter_pvalue_fc(ATACseq_het_vs_WT_BMP_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                           subset_p='padj', subset_fc='log2fc')
ATACseq_hom_vs_het_BMP_2 = filter_pvalue_fc(ATACseq_hom_vs_het_BMP_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                           subset_p='padj', subset_fc='log2fc')

#DORs per comparison in NO
ATACseq_hom_vs_WT_NO_2 = filter_pvalue_fc(ATACseq_hom_vs_WT_NO_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                           subset_p='padj', subset_fc='log2fc')
ATACseq_het_vs_WT_NO_2 = filter_pvalue_fc(ATACseq_het_vs_WT_NO_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                           subset_p='padj', subset_fc='log2fc')
ATACseq_hom_vs_het_NO_2 = filter_pvalue_fc(ATACseq_hom_vs_het_NO_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                           subset_p='padj', subset_fc='log2fc')


#DORs between timepoints
ATACseq_timecomp_BMP_2 = filter_pvalue_fc(ATACseq_timecomp_BMP_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                          subset_p='padj', subset_fc='log2fc')
ATACseq_timecomp_NO_2 = filter_pvalue_fc(ATACseq_timecomp_NO_1, pvalue_threshold=0.01, fc_threshold=1.5,
                                         subset_p='padj', subset_fc='log2fc')

#save the datasets
ATACseq_PAX3_BMP_2.to_csv('ATACseq_PAX3_BMP_DORs.csv')
ATACseq_PAX3_NO_2.to_csv('ATACseq_PAX3_NO_DORs.csv')
ATACseq_timecomp_BMP_2.to_csv('ATACseq_timecomp_BMP_DORs.csv')
ATACseq_timecomp_NO_2.to_csv('ATACseq_timecomp_NO_DORs.csv')

ATACseq_hom_vs_WT_BMP_2.to_csv('ATAseq_hom_vs_WT_BMP.csv')
ATACseq_het_vs_WT_BMP_2.to_csv('ATAseq_het_vs_WT_BMP.csv')
ATACseq_hom_vs_het_BMP_2.to_csv('ATAseq_hom_vs_het_BMP.csv')

ATACseq_hom_vs_WT_NO_2.to_csv('ATAseq_hom_vs_WT_NO.csv')
ATACseq_het_vs_WT_NO_2.to_csv('ATAseq_het_vs_WT_NO.csv')
ATACseq_hom_vs_het_NO_2.to_csv('ATAseq_hom_vs_het_NO.csv')


