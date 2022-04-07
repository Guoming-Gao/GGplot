import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(color_codes=True, style='white')



# batch process all data in a folder
folder1 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Results-20210909-THOR-delta-L941-UGDHOPS'
folder2 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Results-20200315_12dishes'
savingpath = "/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Results-0929colocalization/Supplement_CondensateProperties"
nm_per_pixel = 134
min_tracklength = 5

os.chdir(folder1)
fname_all = [f for f in os.listdir(folder1) if f.endswith('right-condensates.csv')]
lst_df = [pd.read_csv(f) for f in fname_all]

os.chdir(folder2)
fname_all = [f for f in os.listdir(folder2) if f.endswith('right-condensates.csv') and f.find('UGD')>0]
lst_df = lst_df + [pd.read_csv(f) for f in fname_all]

# data without P body, the other HOPS proteins
fname_all = [f for f in os.listdir(folder2) if f.endswith('right-condensates.csv') and f.find('UGD')<0]
lst_df = [pd.read_csv(f) for f in fname_all]

spotInt = np.array([], dtype=float)

df_perfile = lst_df[0]
for df_perfile in lst_df:
    lst_rowID = df_perfile.RowID.unique()

    row = lst_rowID[0]
    for row in lst_rowID:
        df_perrowID = df_perfile[df_perfile.RowID==row]
        lst_tracks = df_perrowID.trackID.unique()

        id = lst_tracks[0]
        for id in lst_tracks:
            df = df_perrowID[df_perrowID.trackID==id]

            if df.shape[0] < min_tracklength:
                continue

            df.keys()
            spotInt = np.append(spotInt, df.totalIntensity.max())




# plotting
plt.figure(figsize=(9, 4), dpi=200)
sns.kdeplot(data=spotInt, fill=True, common_norm=False, color=sns.color_palette()[2])
#sns.kdeplot(data=data, shade=True, log_scale=False, clip=(0.0, 1600))
# plt.title('Radius Distribution (Spot Detector, TrackMate)', fontsize=13, fontweight='bold')
plt.title('Total Intensity Distribution (Spot Detector, TrackMate)', fontsize=13, fontweight='bold')
plt.xlabel('Total Intensity (A.U.)')
plt.legend(["P bodies and HOPS (2x)"])
plt.tight_layout()
fsave = savingpath + '/' + "TotalIntensity_distribution_TrackMate_UGDHOPS.png"
plt.savefig(fsave, format='png')
