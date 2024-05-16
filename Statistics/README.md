# Statistic tests
This file resume the statistical analysis performed by the **stat_test.py** program.
### Statistical analysis
This program allow to test whether samples originate from the same distribution. 
The normality of the distrubutions and equality of variances are assessed. 
Depending on the results, an ANOVA analysis or Kruskal-Wallis test are performed and then post-hoc analysis 
(Dunn's test or Tukey's test).\
The following scheme resume the analysis. 
<img width="852" alt="tests" src="https://github.com/lucasDNS9/Ribes_lab/assets/127426611/3bf51a46-5b7c-467e-bbbf-25d28930b453">

### How to run the program ?
Your data should be a csv file. One column containing the groups you want to compare, 
and another containing the values.\
**Structure :**
<img width="1056" alt="Structure" src="https://github.com/lucasDNS9/Ribes_lab/assets/127426611/4deb61b8-b733-4245-80ef-b11ad3b4e6aa">
The program can be run from the terminal:
```
python3 stat_test.py --data my_data.csv --groupes Groups --values Values --p 0.05
```
# Special Warning
Be careful to the sign you used to delimitate the columns in your dataframe. Need to be ';'
