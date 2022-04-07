import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns
sns.set(color_codes=True, style='white')

###############################################
# Loading Data
fpath1 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210612-dwelltimes.csv'
fpath2 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210909-dwelltimes.csv'
fpath3 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20211007-dwelltimes.csv'
folderpath_save = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD'

df1 = pd.read_csv(fpath1)
df2 = pd.read_csv(fpath2)
df3 = pd.read_csv(fpath3)

# Data cleaning
df1.filename = [name.split('_cell')[0] for name in df1.filename]
df1['RNAtype'] = [name.split('in')[0] for name in df1.filename]
df1['HOPScondition'] = [name.split('in')[1] for name in df1.filename]
df1['date'] = np.repeat('20210612', df1.shape[0])
df1 = df1.replace('miR21ds', 'miR-21ds')
df1 = df1.replace('miR21guidestrand', 'miR-21gs')
df1 = df1.replace('SOX2', "SOX2")
df1 = df1.replace('UGDbefore', '150 mM Na+')
df1 = df1.replace('UGDunderHOPS', '300 mM Na+')

df2.filename = [name.split('-')[1] for name in df2.filename]
df2['RNAtype'] = df2.filename
df2['HOPScondition'] = np.repeat('300 mM Na+', df2.shape[0])
df2['date'] = np.repeat('20210909', df2.shape[0])
df2.RNAtype = df2.RNAtype.replace('dish1','THOR')
df2.RNAtype = df2.RNAtype.replace('dish2','THOR-d')
df2.RNAtype = df2.RNAtype.replace('dish4','L941')

df3.filename = [name.split('_CELL')[0] for name in df3.filename]
df3['RNAtype'] = [name.split(' ')[0] for name in df3.filename]
df3['HOPScondition'] = df3.filename
df3['date'] = np.repeat('20211007', df3.shape[0])
df3.RNAtype = df3.RNAtype.replace('L491','L941')
df3.RNAtype = df3.RNAtype.replace('ACTB','ActB')

for i in range(df3.shape[0]):
    if df3.iloc[i,3][-2] == '1':
        df3.iloc[i,3] = '150 mM Na+'
    elif df3.iloc[i,3][-2] == '2':
        df3.iloc[i,3] = '300 mM Na+'

# df1.to_csv(fpath1[:-4]+'toplot.csv', index=False)
# df2.to_csv(fpath2[:-4]+'toplot.csv', index=False)
# df3.to_csv(fpath3[:-4]+'toplot.csv', index=False)

# Omit Xin's L941 data:
df3 = df3[df3.RNAtype != 'L941']

df_all = pd.concat([df1,df2,df3])
df_all = df_all[df_all.HOPScondition == '300 mM Na+']
df_all.dwelltimes = df_all.dwelltimes.to_numpy(dtype=float)*0.1
df_all = df_all[(df_all.dwelltimes>0) & (df_all.dwelltimes<=20)]

lst_RNA = ['miR-21ds','miR-21gs','THOR','THOR-d','L941','ActB',"SOX2"]

os.chdir(folderpath_save)
i = 0
for i in range(len(lst_RNA)):
    plt.figure(figsize=(9, 3), dpi=200)
    data = df_all[df_all.RNAtype==lst_RNA[i]].dwelltimes
    sns.histplot(data=data, binwidth=0.2, stat='probability', color=sns.color_palette()[i], label=lst_RNA[i])
    plt.legend()
    plt.xlim(0,21)
    plt.ylim(0,0.7)
    plt.xlabel('Dwell Time (s)')

    N = data.size
    comment1 = 'N = ' + str(N)
    plt.text(20.7, 0.52, comment1, horizontalalignment='right', size=12, fontweight='normal')

    percent = round((data[data>=19].size / N) * 100, 2)
    comment2 = 'Static: ' + str(percent) + '%'
    plt.text(20.7, 0.45, comment2, horizontalalignment='right', size=12, fontweight='heavy')

    plt.tight_layout()
    fname_save = 'Dwell Time Distribution-' + lst_RNA[i] + '-staticpercent.png'
    plt.savefig(fname_save, format='png')
    plt.close()
