import os 
import pandas as pd 
import numpy as np 

from tqdm import tqdm
from datetime import date


Path = "Data"
files = os.listdir(Path)

Log2FC = []   # from data
lfcSE = []    # from data 
stat = []     # from data
p_value = []  # from data 
padj = []     # from data 
ensembl = []  # from ?????? 
entrezgene_id = []  # from ?????? 
GeneName = [] # from data
CellLine = [] # from metadata
GSE = []      # from ??????
Rep = []      # from metadata
Dose = []     # from metadata, str
Duration = [] # from metadata, str
minus_log10padj = []  # Compute from padj
Source = []
ids = []

# load meta data
meta_data = pd.read_csv("All_RNAseq_Micorarry_Comparison_metadata_04152022.csv")

index = 0

for file_i in tqdm(files, total=len(files)):
    
    # process each individual file
    mypd = pd.read_csv(os.path.join(Path, file_i))
    mypd.rename(columns={'Unnamed: 0': 'Genename'}, inplace=True)
    assert 'log2FoldChange' in mypd.columns, "log2FoldChange not existed in the dataframe"
    lfcSE_i = 0
    stat_i = 0
    p_value_i = 0
    padj_i = 0
    minus_log10padj_i = 0
    
    # meta data 
    experiment_id_f = ".".join(file_i.split(".")[:-1])
    print(experiment_id_f)
    cellline_f = meta_data[meta_data['Comparison'] == experiment_id_f]['Cell_line'].values[0]
    dose_f = meta_data[meta_data['Comparison'] == experiment_id_f]['Dose'].values[0]
    dose_f = float(dose_f.replace(" ", "").replace("nM", ""))
    duration_f = meta_data[meta_data['Comparison'] == experiment_id_f]['Duration/h'].values[0]
    duration_f = duration_f.replace(",", ";")
    rep_f = meta_data[meta_data['Comparison'] == experiment_id_f]['Replicates'].values[0]
    rep_f = int(rep_f)
    gse_f = meta_data[meta_data['Comparison'] == experiment_id_f]['GSE#'].values[0]
    
    for i in range(len(mypd)):
        ids.append(index)
        index += 1
        data_i = mypd.iloc[i]
        GeneName.append(data_i['Genename'])
        Log2FC.append(data_i['log2FoldChange'])
        if 'lfcSE' in mypd.columns:
            lfcSE_i = data_i['lfcSE']
        if 'stat' in mypd.columns:
            stat_i = data_i['stat']
        if 'p_value' in mypd.columns:
            p_value_i = data_i['p_value']
        if 'padj' in mypd.columns:
            padj_i = data_i['padj']
            if padj_i != 0:
                try:
                    minus_log10padj_i = -np.log10(padj_i)
                except:
                    minus_log10padj_i = 0
        
        lfcSE.append(lfcSE_i)
        stat.append(stat_i)
        p_value.append(p_value_i)
        padj.append(padj_i)
        minus_log10padj.append(minus_log10padj_i)

        # metadata
        CellLine.append(cellline_f)
        Dose.append(dose_f)
        Duration.append(duration_f)
        Rep.append(rep_f)
        Source.append(experiment_id_f)

        # gene id 
        GSE.append(gse_f)



    print(f"Done with {file_i}; Total until now: {index}")




data_pd = pd.DataFrame({
    'id': ids,
    'Log2FC': Log2FC,
    'lfcSE': lfcSE,
    'stat': stat,
    'p_value': p_value,
    'padj': padj,
    'GeneName': GeneName,
    'CellLine': CellLine,
    'GSE': GSE,
    'Rep': Rep,
    'Dose': Dose,
    'Duration': Duration,
    'minus_log10padj': minus_log10padj,
    'Source': Source
})


today = date.today()
d3 = today.strftime("%m%d%y")
data_pd.to_csv(f"RNA_seq_integrate_{d3}.csv", index=False)


    # Log2FC = []   # from data
    # lfcSE = []    # from data 
    # stat = []     # from data
    # p_value = []  # from data 
    # padj = []     # from data 
    # ensembl = []  # from ?????? 
    # entrezgene_id = []  # from ?????? 
    # GeneName = [] # from data
    # CellLine = [] # from metadata
    # GSE = []      # from ??????
    # Rep = []      # from metadata
    # Dose = []     # from metadata, str
    # Duration = [] # from metadata, str
    # minus_log10padj = []  # Compute from padj