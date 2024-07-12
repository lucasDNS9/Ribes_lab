#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 17:51:31 2024

@author: lucasdenis
"""
'''
Section 1: Import
'''
# Import packages
import os
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns

# Define working directory
os.chdir("/Users/lucasdenis/Desktop/ATACseq")

#import dataset
DORs_timecomp_BMP = pd.read_csv("./data/ATACseq_timecomp_BMP_DORs.csv", sep=";")
DORs_timecomp_NO = pd.read_csv("./data/ATACseq_timecomp_NO_DORs.csv", sep=";")

#Define the p-value and fold-change columns
p_comparisons_BMP = ['padj_ctl_d4_vs_ctl_iPS_timecomp', #comparison n°1
                     'padj_ctl_d7_BMP_vs_ctl_d4_timecomp', #comparison n°2
                     'padj_ctl_d7_BMP_vs_ctl_iPS_timecomp'] #comparison n°3

fc_comparisons_BMP = ['log2fc_ctl_d4_vs_ctl_iPS_timecomp', #comparison n°1
                      'log2fc_ctl_d7_BMP_vs_ctl_d4_timecomp', #comparison n°2
                      'log2fc_ctl_d7_BMP_vs_ctl_iPS_timecomp'] #comparaison n°3

p_comparisons_NO = ['padj_ctl_d4_vs_ctl_iPS_timecomp', #comparison n°1
                    'padj_ctl_d7_NO_vs_ctl_d4_timecomp', #comparison n°2
                    'padj_ctl_d7_NO_vs_ctl_iPS_timecomp'] #comparison n°3

fc_comparisons_NO = ['log2fc_ctl_d4_vs_ctl_iPS_timecomp', #comparison n°1
                     'log2fc_ctl_d7_NO_vs_ctl_d4_timecomp', #comparison n°2
                     'log2fc_ctl_d7_NO_vs_ctl_iPS_timecomp'] #comparaison n°3

###############################################################################
'''
Section 2: Classification function
'''
#Define the classification function (with fc_threshold)
def classification(data, p_comparisons, fc_comparisons, 
                      p=0.01, fc=1.5):

    fc_log=math.log2(fc)
    
    profil = ''

    #comparison 1
    if data[p_comparisons[0]] > p:
        profil += '1n'
    else:
        if data[fc_comparisons[0]] > fc_log:
            profil += '1o'
        elif data[fc_comparisons[0]] < -fc_log:
            profil += '1c'
        else:
            profil += '1n'

    #comparison 2
    if data[p_comparisons[1]] > p:
        profil += '2n'
    else:
        if data[fc_comparisons[1]] > fc_log:
            profil += '2o'
        elif data[fc_comparisons[1]] < -fc_log:
            profil += '2c'
        else:
            profil += '2n'

    #comparison 3
    if data[p_comparisons[2]] > p:
        profil += '3n'
    else:
        if data[fc_comparisons[2]] > fc_log:
            profil += '3o'
        elif data[fc_comparisons[2]] < -fc_log:
            profil += '3c'
        else:
            profil += '3n'

    return profil
    

###############################################################################

#Apply the classification function to my dataset BMP
DORs_timecomp_BMP['profils']= DORs_timecomp_BMP.apply(classification, p_comparisons=p_comparisons_BMP,
                                                      fc_comparisons=fc_comparisons_BMP, axis=1)

#Apply the classification function to my dataset NO
DORs_timecomp_NO['profils']= DORs_timecomp_NO.apply(classification, axis=1,
                                                    p_comparisons=p_comparisons_NO,
                                                    fc_comparisons=fc_comparisons_NO)

###############################################################################
'''
Section 3: Analysis and processing
'''
#Count the different elements in the classification
from collections import Counter

# Print the counts of each unique profile
profile_counts_BMP = Counter(list(DORs_timecomp_BMP['profils']))
for profile, count in profile_counts_BMP.items():
    print(f"{profile}: {count}")
    
#Same for NO
profile_counts_NO = Counter(list(DORs_timecomp_NO['profils']))
for profile, count in profile_counts_NO.items():
    print(f"{profile}: {count}")
    
###############################################################################
#Remove the "abberant" profils

junk_profil = ['1o2n3c','1o2o3n','1c2n3o','1c2c3n','1n2c3o']

mask = ~DORs_timecomp_BMP['profils'].isin(junk_profil)
DORs_timecomp_BMP_filt = DORs_timecomp_BMP[mask]

mask_NO = ~DORs_timecomp_NO['profils'].isin(junk_profil)
DORs_timecomp_NO_filt = DORs_timecomp_NO[mask_NO]


###############################################################################
#only keep the average rpkm and profils for the plots
#set index
DORs_timecomp_BMP_filt.set_index('chr:start-end', inplace=True)
DORs_timecomp_NO_filt.set_index('chr:start-end', inplace=True)

#select the columns for subsetting and plotting
rpkm_profils = ['average_RPKM_ctl_iPS','average_RPKM_ctl_d4','average_RPKM_BMP_ctl_d7', 'profils']
DORs_plot_BMP = DORs_timecomp_BMP_filt[rpkm_profils]
rpkm_profils = ['average_RPKM_ctl_iPS','average_RPKM_ctl_d4','average_RPKM_NO_ctl_d7', 'profils']
DORs_plot_NO = DORs_timecomp_NO_filt[rpkm_profils]


###############################################################################
'''
Section 4: Graphs
'''
#Graph function
def plotting_profils(DORs_plot):

    profils = DORs_plot['profils'].unique()
    plots = {}
    
    for profil in profils:
        # Subset the data for the current profile
        subset_DORs = DORs_plot[DORs_plot['profils'] == profil]
        subset_DORs = subset_DORs.drop(columns=['profils'])

        # Calculate means and standard deviations for each timepoint
        medians = subset_DORs.median()
        std_devs = subset_DORs.std()

        timepoints = ['iPSC', 'day 4', 'day 7']

        # Create a DataFrame
        data = pd.DataFrame({
            'Timepoints': timepoints,
            'Medians': medians,
            'Lower': [m - s for m, s in zip(medians, std_devs)],
            'Upper': [m + s for m, s in zip(medians, std_devs)]
        })
        
        # Set the style
        sns.set(style="ticks")

        # Create the plot
        plt.figure(figsize=(5, 4))

        # Plot the means
        sns.lineplot(x='Timepoints', y='Medians', data=data, marker='o')

        # Plot the error bands
        plt.fill_between(data['Timepoints'], data['Lower'], data['Upper'], color='b', alpha=0.2)

        plt.xlabel('')
        plt.ylabel('')
        plt.title(profil)
        
        # Save the plot as a PDF
        pdf_filename = f'plot_{profil}.pdf'
        plt.savefig(pdf_filename)

        # Store the plot in the dictionary
        plots[profil] = plt.gcf()

        # Clear the figure to avoid overlap in the next iteration
        plt.clf()

    return plots


###############################################################################
#Create and save the plots for the 18 profils (with and without the z-score normalization)
plots = plotting_profils(DORs_plot_BMP)
plots = plotting_profils(DORs_plot_NO)








