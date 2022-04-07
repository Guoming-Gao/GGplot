import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns
sns.set(color_codes=True, style='white')
from statannot import add_stat_annotation


fpath1 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210612-All_coloc_percent.csv'
fpath2 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210909-All_coloc_percent.csv'
fpath3 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20211007-All_coloc_percent.csv'

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
df3.HOPScondition = df3.HOPScondition.replace()
for i in range(df3.shape[0]):
    if df3.iloc[i,5][-2] == '1':
        df3.iloc[i,5] = '150 mM Na+'
    elif df3.iloc[i,5][-2] == '2':
        df3.iloc[i,5] = '300 mM Na+'

df1.to_csv(fpath1[:-4]+'toplot.csv', index=False)
df2.to_csv(fpath2[:-4]+'toplot.csv', index=False)
df3.to_csv(fpath3[:-4]+'toplot.csv', index=False)

df_all = pd.concat([df1,df2,df3])

df_PB = deepcopy(df_all[['date','HOPScondition','RNAtype']])
df_PB['CondensateType'] = np.repeat('P bodis', df_PB.shape[0])
df_PB['Percentage of Colocalization'] = df_all.coloc_percent_PB
df_HOPS = deepcopy(df_all[['date','HOPScondition','RNAtype']])
df_HOPS['CondensateType'] = np.repeat('HOPS condensates',df_HOPS.shape[0])
df_HOPS['Percentage of Colocalization'] = df_all.coloc_percent_HOPS
df_all_reshape = pd.concat([df_PB,df_HOPS])
df_all_reshape['RNA'] = df_all_reshape.RNAtype
df_all_reshape['Condensate Type'] = df_all_reshape.CondensateType

data = df_all_reshape[df_all_reshape.HOPScondition == '300 mM Na+']
x = 'RNA'
y = 'Percentage of Colocalization'
hue = 'Condensate Type'

sns.boxplot(data=data, x=x, y=y, hue=hue)
# ax = sns.swarmplot(data=data, x=x, y=y, hue=hue, color='.25')
plt.xticks(rotation=30)
plt.tight_layout()
plt.show()

"""
plt.xlabel(x, fontsize=11)
plt.ylabel('Mass (* multimer state)', fontsize=11)

test_results = add_stat_annotation(ax, data=df_plot, x=x, y=y, order=order,
                                   box_pairs=[("y", "n")],
                                   test='t-test_ind', comparisons_correction=None,
                                   text_format='star', loc='inside', verbose=2)














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
"""
