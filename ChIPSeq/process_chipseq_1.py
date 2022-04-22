import os 
import pandas as pd 
import numpy as np 
from datetime import date

# Get meta data information
meta_data_file = "ChIPseq_Inventory_12012020_YW.csv"
meta_data = pd.read_csv(meta_data_file)

GSE = []
Cellline = []
Duration = [] # char
Dose = [] # float
# fill in from the data bed file
chr_num = []
start = []
end = []
peakid = []
mid = [] # computed by start + end / 2
score = []
Source = []

ids = []

index = 0

file_with_no_reference = []
# Bed file
Folders = ['Bed', 'Bed_02052022', 'Bed_03142022', 'Bed_03152022', 'Bed_04102022']
for folder in Folders:
    files_folder = os.listdir(folder)

    for file_f in files_folder:
        # get the meta data first
        experiment_id_f = file_f.split("_")[0]
        if len(meta_data[meta_data['Experiment_ID'] == experiment_id_f]) == 0:
            file_with_no_reference.append(experiment_id_f)
            print(f"\n{file_f}'s experiment id {experiment_id_f} not found in meta data (folder {folder})\n")
        else:
            gse_f = meta_data[meta_data['Experiment_ID'] == experiment_id_f]['GEO_Accession#'].values[0]
            cellline_f = meta_data[meta_data['Experiment_ID'] == experiment_id_f]['Cell_line'].values[0]
            duration_f = meta_data[meta_data['Experiment_ID'] == experiment_id_f]['Treatment_Duration_h'].values[0]
            duration_f = str(duration_f)

            if duration_f != duration_f:
                duration_f = "0"
            else:
                duration_f = duration_f.replace(" ", "").replace(",", "-")
            dose_f = meta_data[meta_data['Experiment_ID'] == experiment_id_f]['E2_Dose_nM'].values[0]
            if not np.isnan(dose_f):
                dose_f = float(dose_f)
            else:
                dose_f = 0.0
            
            # process data 
            file_data_pd = pd.read_csv(os.path.join(folder, file_f), header=None, sep="\t")
            file_data_pd = file_data_pd.iloc[:, :5]
            file_data_pd.columns  = ['chr', 'start', 'end', 'peakid', 'score']
            for idx_i in range(len(file_data_pd)):
                ids.append(index)
                index += 1
                data_point_i = file_data_pd.iloc[idx_i]
                chr_num.append(data_point_i['chr'])
                start.append(int(data_point_i['start']))
                end.append(int(data_point_i['end']))
                mid.append(int( (data_point_i['end'] + data_point_i['start']) / 2  ))
                
                peakid.append(data_point_i['peakid'])
                score.append(data_point_i['score'])

                # append meta data
                GSE.append(gse_f)
                Cellline.append(cellline_f)
                Duration.append(duration_f)
                Dose.append(dose_f)
                Source.append(experiment_id_f)


            print(f"Done with file {file_f}; total data points {index}")
            # break
print(f"total data points {index}")
data_pd = pd.DataFrame({
    'id': ids,
    'GSE': GSE,
    'Cellline': Cellline,
    'Duration': Duration,
    'Dose': Dose,
    'chr_num': chr_num,
    'start': start,
    'end': end,
    'peakid': peakid,
    'mid': mid,
    'score': score,
    'Source': Source
})


today = date.today()
d3 = today.strftime("%m%d%y")
data_pd.to_csv(f"ChIPSeq_integrate_{d3}.csv", index=False)

print(f"file_with_no_reference : {file_with_no_reference}")
