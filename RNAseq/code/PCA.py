#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 16:13:07 2024

@author: lucasdenis
"""
import os
os.chdir('/Users/lucasdenis/Desktop/RNAseq')

#import datasets
import pandas as pd
import numpy as np
sample_group = pd.read_csv('./datasets/PCA_datasets/sample_group.csv', sep=';')
RNAseq = pd.read_csv('./datasets/PCA_datasets/RNAseq_fpkm_filt.csv', sep=';')
RNAseq_d7 = pd.read_csv('./datasets/PCA_datasets/RNAseq_d7_fpkm_filt.csv', sep=';')

RNAseq_BMP = pd.read_csv('./datasets/PCA_datasets/RNAseq_PAX3_BMP_fpkm_filt.csv', sep=';')
RNAseq_NO = pd.read_csv('./datasets/PCA_datasets/RNAseq_PAX3_NO_fpkm_filt.csv', sep=';')

#set the indexes
RNAseq.set_index('gene_symbol', inplace=True)
RNAseq_d7.set_index('gene_symbol', inplace=True)

RNAseq_BMP.set_index('gene_symbol', inplace=True)
RNAseq_NO.set_index('gene_symbol', inplace=True)

#log transform the data for the PCA (emphasize the most significant sources of variation in the dataset)
RNAseq_log = RNAseq.applymap(np.log1p)
RNAseq_d7_log = RNAseq_d7.applymap(np.log1p)

RNAseq_BMP_log = RNAseq_BMP.applymap(np.log1p)
RNAseq_NO_log = RNAseq_NO.applymap(np.log1p)

###############################################################################
'''
Graphical functions
'''
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

#Scree plot
def scree_plot(data, subset, title=None):
    
    #Select columns that start with subset
    fpkm_data = data[[col for col in data.columns if col.startswith(subset)]]

    # Transpose the DataFrame
    fpkm_data_transposed = fpkm_data.T

    # Perform PCA
    pca = PCA(n_components=int(len(fpkm_data.columns)/2))
    pca.fit_transform(fpkm_data_transposed)

    #percentage of variation
    per_var = np.round(pca.explained_variance_ratio_*100, decimals=1)
    labels = [str(x) for x in range(1,len(per_var)+1)]

    #Plot features
    plt.figure(figsize=(10, 6))
    plt.bar(x=range(1, len(per_var) + 1), height=per_var, tick_label=labels)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Components')
    plt.title('Scree Plot')
    plt.axhline(y=5, color='red', linestyle='--', label='5% Explained Variance')
    plt.legend()
    
    #annotation
    for i in range(3):
        plt.annotate(f'PC{i + 1}: {per_var[i]}%', (i + 1, per_var[i]), 
                     textcoords="offset points", xytext=(0, 10), ha='center')
    
    if title:
        plt.savefig(title, format='pdf')

##########################################


#PCA function
def plot_PCA(data, subset, sample_group, n_components=2, title=None):

    #Select fpkm columns (raw and average)
    fpkm_data = data[[col for col in data.columns if col.startswith(subset)]]
        
    # Transpose the DataFrame
    fpkm_data_transposed = fpkm_data.T
    
    
    ############
    #mapping for graphic representation
    ############

    #mapping the corresponding group to each samples
    key_sample = sample_group['sample_name']
    values_sample = sample_group['group_name']
    sample_mapping = dict(zip(key_sample,values_sample))

    samples = []
    # Iterate over the elements in fpkm_data_transposed.index
    for element in fpkm_data_transposed.index:
        # Check each key in the dictionary
        for key_sample in sample_mapping.keys():
            # If the key is found in the element
            if key_sample in element:
                # Append the corresponding value to the samples list
                samples.append(sample_mapping[key_sample])
                break  # Break the loop once a match is found

    #color mapping
    color_mapping = {
        'ctl_NO' : 'tab:green',
        'ctl_BMP' : 'tab:green',
        'het_NO' : 'tab:orange',
        'het_BMP' : 'tab:orange',
        'hom_NO' : 'tab:red',
        'hom_BMP' : 'tab:red',
        'ctl_iPS' : 'tab:green',
        'ctl_d4' : 'tab:green'}

    # Generate a list of colors corresponding to each condition
    colors = []
    for condition in samples:
        if 'ctl_NO' in condition:
            colors.append(color_mapping['ctl_NO'])
        elif 'ctl_BMP' in condition:
            colors.append(color_mapping['ctl_BMP'])
        elif 'het_NO' in condition:
            colors.append(color_mapping['het_NO'])
        elif 'het_BMP' in condition:
            colors.append(color_mapping['het_BMP'])
        elif 'hom_NO' in condition:
            colors.append(color_mapping['hom_NO'])
        elif 'hom_BMP' in condition:
            colors.append(color_mapping['hom_BMP'])
        elif 'ctl_iPS' in condition:
            colors.append(color_mapping['ctl_iPS'])        
        elif 'ctl_d4' in condition:
            colors.append(color_mapping['ctl_d4'])      
            
    #shape_mapping
    shape_mapping = {
        'BMP': 'o',
        'NO': 'x',
        'd4': 'd',
        'iPS': '*'
    }

    shapes = []
    for condition in samples:
        if 'BMP' in condition:
            shapes.append(shape_mapping['BMP'])
        elif 'NO' in condition:
            shapes.append(shape_mapping['NO'])
        elif 'd4' in condition:
            shapes.append(shape_mapping['d4'])
        elif 'iPS' in condition:
            shapes.append(shape_mapping['iPS'])
   
    
    # Perform PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(fpkm_data_transposed)
    
    #percentage of variation
    per_var = np.round(pca.explained_variance_ratio_*100, decimals=1)

    # Create a DataFrame for the PCA results
    pca_RNAseq = pd.DataFrame(data=pca_result, 
                              columns=[f'PC{i+1}' for i in range(n_components)])

    # Set the index to the list of samples
    pca_RNAseq.set_index(pd.Index(fpkm_data_transposed.index), inplace=True)

    # Plot the PCA results
    plt.figure(figsize=(8, 6))
    
    # Cosmetic PCA plot
    for i, (pc1, pc2) in enumerate(zip(pca_RNAseq['PC1'], pca_RNAseq['PC2'])):
        plt.scatter(pc1, pc2, color=colors[i], marker=shapes[i], s=60, label=sample_group['group_name'].iloc[i])
    
    plt.title('', fontsize=14)
    
    if n_components == 2:
        plt.xlabel('PC 1 ('+str(per_var[0])+'%)', fontsize=14)
        plt.ylabel('PC 2 ('+str(per_var[1])+'%)', fontsize=14)

    else:
        print(f"PCA Results for {n_components} components:\n", pca_RNAseq)

    if title:
        plt.savefig(title, format='pdf')


###############################################################################
'''
Graphs
'''

#Scree plot
scree_plot(RNAseq_log, subset='fpkm')
scree_plot(RNAseq_d7_log, subset='fpkm')

#Perform the PCA (on fpkm Z-score normalized)
plot_PCA(data=RNAseq_log, subset='fpkm', sample_group=sample_group, n_components=2, title=None)
plot_PCA(data=RNAseq_d7_log, subset='fpkm', sample_group=sample_group, n_components=2, title=None)

plot_PCA(data=RNAseq_BMP, subset='fpkm', sample_group=sample_group, n_components=2, title=None)
plot_PCA(data=RNAseq_NO, subset='fpkm', sample_group=sample_group, n_components=2, title=None)



