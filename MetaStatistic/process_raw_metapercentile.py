import pandas as pd 
import os 


file_pd = pd.read_csv("Raw_RNA_Micro_statistics.csv")
file_pd.columns = ["Gene"] + list(file_pd.columns[1:])

# process this into tidy format
file_pd_melt = pd.melt(file_pd, id_vars="Gene", var_name="Source", value_name="Percentile")
file_pd_melt = file_pd_melt.dropna().sort_values('Gene')
file_pd_melt['id'] = file_pd_melt.index
print(file_pd_melt)
file_pd_melt.to_csv("gene_raw_meta_percentile.csv", index=False)

