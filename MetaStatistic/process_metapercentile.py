import os 
import pandas as pd 
import numpy as np 

file_pd = pd.read_csv("gene_meta_percentile.csv")
file_pd.columns = ['Gene', 'bin', 'count']

file_pd.to_csv("gene_meta_percentile.csv")