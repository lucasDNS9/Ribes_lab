#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 17:30:24 2024

@author: lucasdenis
"""
###############################################################################

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from scipy.stats import kruskal, levene, f_oneway, shapiro
from scikit_posthocs import posthoc_dunn
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import sys
from io import StringIO
import argparse
import pandas as pd

###############################################################################
#function to generate a pdf :
def generate_pdf(output_text, filename):
    
    c = canvas.Canvas(filename, pagesize=letter)
    text = c.beginText(100, 750)
    text.setFont("Helvetica", 12)
    for line in output_text.split('\n'):
        text.textLine(line)
    c.drawText(text)
    c.save()


###############################################################################
#Function to perform the statistical tests
def test_stat(data, groupes, values, p=0.05):
    
    #allow to generate a pdf with the stat results
    output = StringIO()
    sys.stdout = output

    #remove the row with NaN
    filtered_data = data[data[values].notna()]
    
    grouped_data = data.groupby(groupes)[values].apply(list)
    
    #Vérification de la normalité des distributions : Shapiro-test
    p_values_Shapiro_test=[]

    print('Shapiro test :')
    for x, y in grouped_data.items():
        # Perform Shapiro-Wilk test on the values within each group
        stat_Shapiro, p_Shapiro = shapiro(y)
        print(f"Genotype: {x}")
        print("Test statistic:", stat_Shapiro)
        print("P-value:", p_Shapiro)
        
        p_values_Shapiro_test.append(p_Shapiro)
       
    if any(x < p for x in p_values_Shapiro_test): #do the Kruskal test
        
        print('Conclusion : The data are not normaly distributed.\n')
        print('We perform a Kruskal-Wallis test :')
    
        stat_Kruskal, p_Kruskal = kruskal(*grouped_data)
        # Afficher les résultats
        print("Statistique de test de Kruskal-Wallis :", stat_Kruskal)
        print("Valeur de p :", p_Kruskal)
        
        if p_Kruskal < p: 
            
            print('Conclusion : There is a difference between groups.\n')
            print("We perform the post-hoc Dunn's test :")
            
            dunn_results = posthoc_dunn(filtered_data, val_col=values, group_col=groupes)
            print("Dunn's Test results :")
            print(dunn_results)
        
        else :
            print('Conclusion : There is no statistical differences between groups.\n')

    else : #test the homogeneity of the variance across group
        print('Conclusion : The data are normaly distributed.\n')
        print('We perform the Levene test to assess the homogeneity of variances :')
        
        stat_Levene, p_Levene = levene(*grouped_data)
        # Print the test results
        print("Levene Test Statistic:", stat_Levene)
        print("p-value:", p_Levene)

        # Interpret the results
        if p_Levene < p:
            print("Conclusion : Variances are not equal.\n")
            print('We perform a Kruskal-Wallis test :')
        
            stat_Kruskal, p_Kruskal = kruskal(*grouped_data)
            # Afficher les résultats
            print("Statistique de test de Kruskal-Wallis :", stat_Kruskal)
            print("Valeur de p :", p_Kruskal)
            
            if p_Kruskal < p: 
                
                print('Conclusion : There is a differences between groups.\n')
                print("We perform the post-hoc Dunn's test :")
                
                dunn_results = posthoc_dunn(filtered_data, val_col=values, group_col=groupes)
                print("Dunn's Test results :")
                print(dunn_results)
            
            else :
                print('Conclusion : There is no statistical differences between groups.\n')

        else:
            print("Conclusion : Variances are equal.\n")
            print('We perform an ANOVA test :')
            
            # Perform ANOVA
            stat_ANOVA, p_ANOVA = f_oneway(*grouped_data)
            # Print the test results
            print("ANOVA F-statistic:", stat_ANOVA)
            print("p-value:", p_ANOVA)

            if p_ANOVA < p:
                print("Conclusion: There is a significant difference between group means.\n ")
                print('We perform a Tukey post-hoc test :')
                
                # Perform Tukey post-hoc test
                tukey_results = pairwise_tukeyhsd(filtered_data[values], filtered_data[groupes])
                print("Tukey HSD results :")
                print(tukey_results)
            
            else:
                print("Conclusion : There is no significant difference between group means.\n")
    
    group_counts = filtered_data[groupes].value_counts()
    # Display unique groups and their counts
    print('number of n:')
    print(group_counts)

    #Get the printed output as a string
    output_text = output.getvalue()
    # Reset stdout
    sys.stdout = sys.__stdout__

    #Generate a pdf file with the stat results
    filename = 'Stat_'+values+'.pdf'
    generate_pdf(output_text, filename)


###############################################################################
#Parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', help='Path to the CSV file containing the data')
    parser.add_argument('--groupes', help='Name of the column containing group labels')
    parser.add_argument('--values', help='Name of the column containing values to compare')
    parser.add_argument('--p', default=0.05, type=float, help='Significance level')
    args = parser.parse_args()
    
    # Load data from CSV file
    data = pd.read_csv(args.data, sep=';')
    
    # Call your function with the provided arguments
    test_stat(data, args.groupes, args.values, args.p)


###############################################################################
#python3 stat_test.py --data data_quantif_CM34_HUCD_SOX2_r.csv --groupes Genotype --values ratio_HUCD --p 0.05







