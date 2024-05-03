#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:59:15 2024

@author: lucasdenis
"""
import os
os.chdir('/Users/lucasdenis/Desktop/RNAseq')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#import datasets
RNAseq_PAX3_BMP_filt_fpkm = pd.read_csv('./datasets/fpkm_filtered/RNAseq_PAX3_BMP_fpkm_filt.csv')
RNAseq_PAX3_NO_filt_fpkm = pd.read_csv('./datasets/fpkm_filtered/RNAseq_PAX3_NO_fpkm_filt.csv')

#set the indexes
RNAseq_PAX3_BMP_filt_fpkm.set_index('gene_symbol', inplace=True)
RNAseq_PAX3_NO_filt_fpkm.set_index('gene_symbol', inplace=True)

###############################################################################


'''
Define the functions
'''
###############################################################################
#function to select one comparison for the volcano plot
def select_comp(data, condition_1, condition_2):
    
    selected_data = pd.concat([data.filter(like='log2fc'),
                           data.filter(like='padj')],axis=1)
    selected_data = selected_data.filter(like=condition_1).filter(like=condition_2)
    selected_data.columns.values[0] = 'log2fc'
    selected_data.columns.values[1] = 'padj'
    
    return selected_data


###############################################################################
#function to extract DEGs from a single comparison
def gene_extraction(data, condition_1, condition_2, fc_threshold, 
                    pvalue_threshold, save=False):
    
    selected_data = select_comp(data, condition_1, condition_2)
    
    fc_threshold_log = np.log2(fc_threshold)
    
    #Select the genes from the new dataset
    mask = (selected_data['padj'] < pvalue_threshold) & (np.abs(selected_data['log2fc']) > fc_threshold_log)
    diff_expr_genes = selected_data[['log2fc', 'padj']][mask]

    if save:
        diff_expr_genes.to_csv('DEGs_'+str(condition_1)+'_vs_'+str(condition_2)+'.csv', 
                               index=False)
    
    return diff_expr_genes


###############################################################################
#Volcano plot function
def volcano_plot(data, condition_1, condition_2, fc_threshold, 
                 pvalue_threshold=0.01, label=False, save=False, pdf=False):
    
    #select data for the comparison (with the select_comp function)
    comp = select_comp(data, condition_1, condition_2)
    
    #transform the pvalue
    comp['padj'] = comp['padj'].apply(lambda x: -np.log10(x))
    
    # Thresholds for significance and fold change
    pvalue_lim_log = -np.log10(pvalue_threshold)
    fc_lim_log2 = np.log2(fc_threshold)
    
    #Create a volcano plot
    plt.figure(figsize=(10, 6))
    plt.scatter(comp['log2fc'], comp['padj'], color='grey', 
                alpha=0.7, edgecolors='none')
    
    #Set up the axis
    plt.xlim(-11,5)
    plt.ylim(-1,25)

    # Highlight significant points
    mask_up = (comp['padj'] > pvalue_lim_log) & (comp['log2fc'] > fc_lim_log2)
    plt.scatter(comp['log2fc'][mask_up], comp['padj'][mask_up], 
                color=(0.949217, 0.517763, 0.295662, 1.0), alpha=0.5)
    plt.text(3, 22, 'n='+str(sum(mask_up))+' genes', fontsize=16, 
             ha='center', color=(0.949217, 0.517763, 0.295662, 1.0))
    
    mask_down = (comp['padj'] > pvalue_lim_log) & (comp['log2fc'] < -fc_lim_log2)
    plt.scatter(comp['log2fc'][mask_down], comp['padj'][mask_down], 
                color=(0.562738, 0.051545, 0.641509, 1.0), alpha=0.5)
    plt.text(-9, 22, 'n='+str(sum(mask_down))+' genes', fontsize=16, 
             ha='center', color=(0.562738, 0.051545, 0.641509, 1.0))


    if save:
        #Save a pdf with the gene differentially expressed
        genes_diff = gene_extraction(data, condition_1, condition_2, fc_threshold, pvalue_threshold)
        genes_diff.to_excel(str(condition_1)+'vs_'+str(condition_2)+'.xlsx', index=True)
    
    # Set labels and title
    plt.xlabel('log2(fold change)')
    plt.ylabel('-log(p-value)')
    plt.title(str(condition_1)+' vs. '+ str(condition_2))
    plt.axhline(y=pvalue_lim_log, color='black', linestyle='--', linewidth=1, 
                label='Significance Threshold: '+str((pvalue_threshold)))
    plt.axvline(x=fc_lim_log2, color='blue', linestyle='--', linewidth=1, 
                label='Fold Change Threshold: '+str((fc_threshold)))
    plt.axvline(x=-fc_lim_log2, color='blue', linestyle='--', linewidth=1)

    # Add legend
    #plt.legend()
    
    if pdf:
        plt.savefig(str(condition_1)+'_vs_'+str(condition_2)+'.pdf', format='pdf')
        
###############################################################################
'''
Volcano generation (important to give data already filtered for fpkm values)
'''
#BMP hom vs ctl
volcano_plot(RNAseq_PAX3_BMP_filt_fpkm, condition_1='hom_BMP_D7', 
             condition_2='ctl_BMP_D7', fc_threshold=1.5, pvalue_threshold=0.01,
             label=False, pdf=True, save=True)

#BMP het vs ctl
volcano_plot(RNAseq_PAX3_BMP_filt_fpkm, condition_1='het_BMP_D7', 
             condition_2='ctl_BMP_D7', fc_threshold=1.5, pvalue_threshold=0.01,
             label=False, pdf=True, save=True)

#BMP hom vs het
volcano_plot(RNAseq_PAX3_BMP_filt_fpkm, condition_1='hom_BMP_D7', 
             condition_2='het_BMP_D7', fc_threshold=1.5, pvalue_threshold=0.01,
             label=False, pdf=True, save=True)


#NO BMP hom vs ctl
volcano_plot(RNAseq_PAX3_NO_filt_fpkm, condition_1='hom_NO_D7', 
             condition_2='ctl_NO_D7', fc_threshold=1.5, pvalue_threshold=0.01,
             label=False, pdf=True, save=True)

#NO BMP het vs ctl
volcano_plot(RNAseq_PAX3_NO_filt_fpkm, condition_1='het_NO_D7', 
             condition_2='ctl_NO_D7', fc_threshold=1.5, pvalue_threshold=0.01,
             label=False, pdf=True, save=True)
    
#NO BMP hom vs het
volcano_plot(RNAseq_PAX3_NO_filt_fpkm, condition_1='hom_NO_D7', 
             condition_2='het_NO_D7', fc_threshold=1.5, pvalue_threshold=0.01,
             label=False, pdf=True, save=True)





        