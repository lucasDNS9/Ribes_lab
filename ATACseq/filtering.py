#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 10:34:25 2024

@author: lucasdenis
"""
import numpy as np

"""
Filtering functions :
"""
#function to filter on fpkm value
def filter_rpkm(data, rpkm_threshold, subset='average_RPKM_'):

    # Identify columns containing "average_fpkm" in their names
    columns_to_check = [col for col in data.columns if col.startswith(subset)]
    
    #filter rows
    filtered_data = data[data[columns_to_check].gt(rpkm_threshold).any(axis=1)]

    return filtered_data


#function to filter on pvalue
def filter_pvalue(data, pvalue_threshold, subset='padj'):

    # Identify columns containing subset in their names
    columns_to_check = data.filter(like=subset).columns
    
    #filter rows
    filtered_data = data[data[columns_to_check].lt(pvalue_threshold).any(axis=1)]

    return filtered_data


#function to filter on fold-change
def filter_fc(data, fc_threshold, subset='log2fc'):

    # Identify columns containing "log2fc" in their names
    columns_to_check = data.filter(like=subset).columns
    
    #filter rows > threshold or < -threshold
    filtered_data = data[(data[columns_to_check] > np.log2(fc_threshold)).any(axis=1)| 
                         (data[columns_to_check] < np.log2(1/fc_threshold)).any(axis=1)]

    return filtered_data


#filter rows on pvalue and fc
def filter_pvalue_fc(data, pvalue_threshold, fc_threshold, subset_p='padj', subset_fc='log2fc'):
    # Identify columns containing subset in their names
    columns_to_check_p = data.filter(like=subset_p).columns
    columns_to_check_fc = data.filter(like=subset_fc).columns
    
    filtered_data_list = []  # Initialize an empty list to store filtered data
    
    # Iterate over each comparison
    for i in range(len(columns_to_check_p)):
        pvalue_column = columns_to_check_p[i]
        fc_column = columns_to_check_fc[i]
        
        #  Extract DEGs for every comparisons
        filtered_data_i = data[(data[pvalue_column] < pvalue_threshold) & 
                               ((data[fc_column] > np.log2(fc_threshold)) | 
                                (data[fc_column] < np.log2(1/fc_threshold)))].index.tolist()
        
        # create a list with all the DEGs
        filtered_data_list.extend(filtered_data_i)
    
    # remove duplicate DEGs
    unique_indices = list(set(filtered_data_list))
    
    #Create a mask from DEGs list
    mask = data.index.isin(unique_indices)

    # Select rows based on the mask
    filtered_data = data[mask]
    
    return filtered_data


