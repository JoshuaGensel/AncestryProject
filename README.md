# AncestryProject

Collaboratory Popgen project of Vladimir Bajic and Leon Joshua Gensel

A detailed report of this project can be found [here](https://rpubs.com/JGen/AncestryProject).

## Used Programs + Versions:

SLiM = 3.7.1 (conda installation)

python = 3.8.5,
tskit = 0.4.1,
biopython = 1.79,
ete3 = 3.1.2,
msprime = 1.1.0,
pyslim = 0.700,
numpy = 1.19.2,

R = 4.2.1,
tidyverse = 1.3.2,
ggplot2 = 3.3.6,
reshape2 = 1.4.4,
betareg = 3.1-4,
imager = 0.42.13,

## How to run the program

To run the entire programm you just need to clone the repository and run the 'automation.py' file 
from the command line (given you have the used software from above installed). E.g. you could clone 
the repository to your local machine, open a conda shellin the src directory and execute 
'python automation.py -n 10' to start 10 consecutive runs of the Simulations + processing pipeline. 
This will append to the treedata .csv-files in '/data/tree_analysis_data'. If you only want your 
own data, delete the .csv files before hand. If you want to change any parameters you can edit the 
source code in automation.py, except for changes in TA. These have to be edited in the 
'SLiM_model.slim' script. The report.Rmd is a reproducabile R Markdown report. In case you ran 
your own simulations the report can still be compiled and will use your data.