#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 21:54:09 2024

@author: lucasdenis
"""
'''
Section 1: Import
'''
#import packages
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#Define working directory
os.chdir("/Users/lucasdenis/Desktop/ATACseq")

#import DORs from time comparisons and the associated profil
DORs_timecomp = pd.read_csv('./results/DORs_timecomp_profils/DORs_timecomp_profils.csv', sep=';')
DORs_timecomp.set_index('chr:start-end', inplace=True)

#import DORs hom_vs_WT_BMP and associated genes
DORs_data = pd.read_csv('./results/lists_DORs-genes_association/ATACseq_hom_vs_WT_BMP_DORs_genes.csv', sep=';')

#import DEGs hom_vs_WT_BMP from RNAseq analysis
DEGs_data = pd.read_csv('./results/from_RNAseq/DEGs_hom_BMP_D7_vs_ctl_BMP_D7.csv', sep=';')

#cross DEGs list with the gene associated with DORs
DEGs_DORs_genes = DORs_data[DORs_data['gene'].isin(DEGs_data['gene_symbol'])]
DEGs_DORs_genes.to_csv('data.csv')
list_DEGs_DORs_genes = list(DEGs_DORs_genes['chrom:start-end'])


###############################################################################
'''
Section 2: plot function
'''
#Function to draw the graph and higlight regions
def plotting_profils_list(DORs_plot, list_regions):

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

        # Count the number of regions to be highlighted
        highlighted_regions = 0

        # Highlight regions if present in DORs_PAX3_loc
        for region in subset_DORs.index:
            if region in list_regions:
                # Plot the line for this region
                plt.plot(timepoints, subset_DORs.loc[region], color='red', linewidth=1)
                highlighted_regions += 1
                # Print the details of the highlighted regions
                print(f'Highlighted region: {region}, Profile: {profil}')
        
        # Print the total number of highlighted regions for this profile
        print(f'number of highlighted regions for profile {profil}: {highlighted_regions}')
        
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

################
#apply the function
plot_DORs_DEGs = plotting_profils_list(DORs_plot=DORs_timecomp, list_regions=list_DEGs_DORs_genes)



