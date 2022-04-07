import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
from scipy.optimize import curve_fit
import seaborn as sns
sns.set(color_codes=True, style='white')

###############################################
# Loading Data
fpath_double = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Dwell Time fitting results_double.csv'
fpath_single = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Dwell Time fitting results_single.csv'
folderpath_save = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD'

df_double = pd.read_csv(fpath_double)
df_single = pd.read_csv(fpath_single)

lst_RNA = ['miR-21ds','miR-21gs','THOR','THOR-d','L941','ActB',"SOX2"]

os.chdir(folderpath_save)
plt.figure(figsize=(5, 5), dpi=200)
i=0
for i in range(len(lst_RNA)):
    row = df_single[df_single.RNA==lst_RNA[i]]
    if row.R2.to_numpy(dtype=float)[0] > 0.95:
        plt.scatter(100, row.halflife, color=sns.color_palette()[i], s=50, label=lst_RNA[i])
    else:
        row = df_double[df_double.RNA==lst_RNA[i]]
        plt.scatter([row.percent1,row.percent2],[row.halflife1,row.halflife2], color=sns.color_palette()[i],s=50,label=lst_RNA[i])

plt.legend()
plt.xlabel('Percentage (%)')
plt.ylabel('Dwelling Half-Life (s)')
plt.tight_layout()
fname_save = 'Dwelling Half-Life distribution among RNAs.png'
plt.savefig(fname_save, format='png')
plt.close()
