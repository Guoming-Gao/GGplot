import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns
sns.set(color_codes=True, style='white')
from statannot import add_stat_annotation


fpath1 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210612-All_coloc_percent_percellroi.csv'
fpath2 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210909-All_coloc_percent_percellroi.csv'
fpath3 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20211007-All_coloc_percent_percellroi.csv'
folderpath_save = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD'

df1 = pd.read_csv(fpath1)
df2 = pd.read_csv(fpath2)
df3 = pd.read_csv(fpath3)

df1.filename = [name.split('_cell')[0] for name in df1.filename]
df1['RNAtype'] = [name.split('in')[0] for name in df1.filename]
df1['HOPScondition'] = [name.split('in')[1] for name in df1.filename]
df1['date'] = np.repeat('20210612', df1.shape[0])
df1 = df1.replace('miR21ds', 'miR-21 double strand')
df1 = df1.replace('miR21guidestrand', 'miR-21 guide strand')
df1 = df1.replace('SOX2', 'SOX2 mRNA')
df1 = df1.replace('UGDbefore', '150 mM Na+')
df1 = df1.replace('UGDunderHOPS', '300 mM Na+')

df2.filename = [name.split('-')[1] for name in df2.filename]
df2['RNAtype'] = df2.filename
df2['HOPScondition'] = np.repeat('300 mM Na+', df2.shape[0])
df2['date'] = np.repeat('20210909', df2.shape[0])
df2.RNAtype = df2.RNAtype.replace('dish1','THOR lncRNA')
df2.RNAtype = df2.RNAtype.replace('dish2','THOR-delta lncRNA')
df2.RNAtype = df2.RNAtype.replace('dish4','L941 lncRNA')

df3.filename = [name.split('_CELL')[0] for name in df3.filename]
df3['RNAtype'] = [name.split(' ')[0] for name in df3.filename]
df3['HOPScondition'] = df3.filename
df3['date'] = np.repeat('20211007', df3.shape[0])
df3.RNAtype = df3.RNAtype.replace('L491','L941 lncRNA')
df3.RNAtype = df3.RNAtype.replace('ACTB','beta-Actin mRNA')

for i in range(df3.shape[0]):
    if df3.iloc[i,6][-2] == '1':
        df3.iloc[i,6] = '150 mM Na+'
    elif df3.iloc[i,6][-2] == '2':
        df3.iloc[i,6] = '300 mM Na+'

df1.to_csv(fpath1[:-4]+'toplot.csv', index=False)
df2.to_csv(fpath2[:-4]+'toplot.csv', index=False)
df3.to_csv(fpath3[:-4]+'toplot.csv', index=False)

df_all = pd.concat([df1,df2,df3])
df_all['RNA'] = df_all.RNAtype
df_all['Percentage of Colocalization'] = df_all.percent
df_all['Number of RNAs per cell'] = df_all['N_total in roi']


data = df_all[df_all.HOPScondition == '300 mM Na+']
x = 'RNA'
y = 'Percentage of Colocalization'
order = ['miR-21 double strand', 'miR-21 guide strand',
       'THOR lncRNA', 'THOR-delta lncRNA', 'L941 lncRNA',
       'beta-Actin mRNA', 'SOX2 mRNA']
box_pairs = [
    ('miR-21 double strand', 'miR-21 guide strand'),
    ('THOR lncRNA', 'THOR-delta lncRNA'),
    ('THOR lncRNA', 'L941 lncRNA'),
    ('beta-Actin mRNA', 'SOX2 mRNA')
]

plt.figure(figsize=(5, 6), dpi=200)
ax = sns.boxplot(data=data, x=x, y=y, order=order)
ax = sns.swarmplot(data=data, x=x, y=y, order=order, color='.25')
test_results = add_stat_annotation(ax, data=data, x=x, y=y, order=order,
                                   box_pairs=box_pairs,
                                   test='t-test_ind', comparisons_correction=None,
                                   text_format='star', loc='inside', verbose=2)
plt.xticks(rotation=30)
plt.title('RNA Colocalization with HOPS Condensates (per cell)')
plt.tight_layout()
os.chdir(folderpath_save)
plt.savefig('RNA Colocalization with HOPS Condensates-percell.png', format='png')
plt.close()


y = 'Number of RNAs per cell'
plt.figure(figsize=(5, 6), dpi=200)
sns.boxplot(data=data, x=x, y=y, order=order, showfliers = False)
plt.xticks(rotation=30)
plt.title('Number of RNAs per cell')
plt.tight_layout()
os.chdir(folderpath_save)
plt.savefig('Number of RNAs-percell.png', format='png')
