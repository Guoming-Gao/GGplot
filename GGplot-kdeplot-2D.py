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


i = 0
for i in range(len(fname_all)):
    df = pd.read_csv(fname_all[i])
    df = df[df.closest.astype(float)<0.5]

    df_plot1 = pd.DataFrame()
    df_plot2 = pd.DataFrame()
    df_plot1['Condensate Class'] = df['Condensate Class']
    df_plot2['Condensate Class'] = df['Condensate Class']
    df_plot1['log10closest'] = df.closest.astype(float)
    df_plot2['log10closest'] = df.closest.astype(float)
    df_plot1['dwell'] = df.dwell_time_s.astype(float)
    # df_plot2['log10D'] = np.log10(df.dropna().D_um2_s.astype(float))
    df_plot2['log10D'] = np.log10(df.D_um2_s.astype(float))


    title = lst_title[i]


    plt.figure(dpi=200)
    g = sns.kdeplot(data=df_plot1, x='log10closest', y='dwell', hue='Condensate Class', alpha=.7, fill=True, levels=[0.001, 0.005, 0.01, 0.1, 0.3, 0.7, 1])
    g.get_legend().set_title(None)
    g.get_legend().set_bbox_to_anchor((0.22,1))
    plt.title(title)
    plt.xlabel('log10(closest distance, normalized)')
    plt.ylabel('Dwell time (s)')
    plt.tight_layout()
    plt.savefig(title+'-closest-dwell.png', format='png')

    plt.figure(dpi=200)
    g = sns.kdeplot(data=df_plot2, x='log10closest', y='log10D', hue='Condensate Class', alpha=.7, fill=True, levels=[0.001, 0.005, 0.01, 0.1, 0.3, 0.7, 1])
    g.get_legend().set_title(None)
    g.get_legend().set_bbox_to_anchor((0.22,1))
    plt.title(title)
    plt.xlabel('log10(closest distance, normalized)')
    plt.ylabel('log10(D, um^2/s)')
    plt.tight_layout()
    plt.savefig(title+'-closest-log10D.png', format='png')

lst = df_plot1['Condensate Class'].unique()
plt.figure(dpi=200)
ax = plt.gca()
for j in list(range(2)):
    label = lst[j]
    print(label)
    data = df_plot2[df_plot2['Condensate Class']==label]
    sns.kdeplot(data=df_plot2, x='log10closest', y='log10D', alpha=.7, fill=True, color=sns.color_palette()[j], label=label, ax=ax)
    #plt.xlim(-5,1)
    #plt.ylim(-5,1)
plt.title(title)
plt.xlabel('log10(closest distance, normalized)')
plt.ylabel('log10(D, um^2/s)')
plt.legend()
plt.tight_layout()
plt.savefig(title+'-closest-log10D.png', format='png')
