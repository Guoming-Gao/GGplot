import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
sns.set(color_codes=True, style='white')

# batch process all data in a folder
folderpath = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Results-0929colocalization'

# main:
os.chdir(folderpath)
# fname_all = [f.split('-') for f in os.listdir(folderpath) if f.endswith('-Results.csv')]
# fname_start = list(set(['-'.join(map(str, f[:-2])) for f in fname_all]))

fname_all = [f for f in os.listdir(folderpath) if f.endswith('-Results.csv') and f.find('UGD')>0]
lst_title = [f.split('-Results')[0] for f in os.listdir(folderpath) if f.endswith('-Results.csv') and f.find('UGD')>0]
# fname_start = ['-'.join(f.split('-')[:-1]) for f in fname_all]

percent_P = np.array([], dtype=float)
percent_HOPS = np.array([], dtype=float)

for i in range(len(fname_all)):
    df = pd.read_csv(fname_all[i])
    df['Condensate Class'].unique()
    N_all = df.shape[0]
    df = df[df.dwell_time_s.astype(float)>0.3]


    N_P = df[df['Condensate Class']=='P body'].shape[0]
    N_HOPS = df[df['Condensate Class']=='HOPS'].shape[0]
    percent_P = np.append(percent_P, N_P / N_all * 100)
    percent_HOPS = np.append(percent_HOPS, N_HOPS / N_all * 100)

lst_title
percent_P
percent_HOPS.size
ratio_HOPStoP = percent_HOPS/percent_P

selector = [f.startswith('PP7') for f in lst_title]

df1 = pd.DataFrame()
df1['title'] = np.append(np.array(lst_title, dtype=object)[selector], np.array(lst_title, dtype=object)[selector])
df1['percent'] = np.append(percent_HOPS[selector], percent_P[selector])
df1['Condensate Class'] = np.append(np.repeat('HOPS',3), np.repeat('P body',3))
plt.figure(figsize=(4, 6), dpi=200)
sns.barplot(x="title", y="percent", hue="Condensate Class", data=df1)
plt.ylabel('Interaction Percentage among All RNAs (%)')
plt.xticks(rotation=45)
plt.xlabel('')
plt.tight_layout()
plt.savefig('barplot-percent-PP7.png', format='png')

df2 = pd.DataFrame()
df2['title'] = np.array(lst_title, dtype=object)[selector]
df2['ratio'] = ratio_HOPStoP[selector]
plt.figure(figsize=(4, 6), dpi=100)
sns.barplot(x="title", y="ratio", data=df2)
plt.ylabel('RNA-HOPS to RNA-P body Interaction Ratio')
plt.xticks(rotation=45)
plt.xlabel('')
plt.tight_layout()
plt.savefig('barplot-ratio-PP7.png', format='png')



selector = [f.startswith('beadloading') for f in lst_title]

df1 = pd.DataFrame()
df1['title'] = np.append(np.array(lst_title, dtype=object)[selector], np.array(lst_title, dtype=object)[selector])
df1['percent'] = np.append(percent_HOPS[selector], percent_P[selector])
df1['Condensate Class'] = np.append(np.repeat('HOPS',3), np.repeat('P body',3))
plt.figure(figsize=(4, 6), dpi=200)
sns.barplot(x="title", y="percent", hue="Condensate Class", data=df1)
plt.ylabel('Interaction Percentage among All RNAs (%)')
plt.xticks(rotation=45)
plt.xlabel('')
plt.tight_layout()
plt.savefig('barplot-percent-beadloading.png', format='png')

df2 = pd.DataFrame()
df2['title'] = np.array(lst_title, dtype=object)[selector]
df2['ratio'] = ratio_HOPStoP[selector]
plt.figure(figsize=(4, 6), dpi=100)
sns.barplot(x="title", y="ratio", data=df2)
plt.ylabel('RNA-HOPS to RNA-P body Interaction Ratio')
plt.xticks(rotation=45)
plt.xlabel('')
plt.tight_layout()
plt.savefig('barplot-ratio-beadloading.png', format='png')
