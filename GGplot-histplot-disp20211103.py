import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns
sns.set(color_codes=True, style='white')

###############################################
# Loading Data
folderpath = '/Volumes/AnalysisGG/KNIME_results_from_greatlakes/20211007-RNAs'
folderpath_save = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD'
os.chdir(folderpath)
lst_files = [f for f in os.listdir(folderpath) if f.endswith("_totaldisp.csv")]
common_names = np.unique([f.split('_CELL')[0] for f in lst_files])

lst_df_1x = [pd.read_csv(f) for f in lst_files if f.startswith('ACTB  BEFORE-1X')]
lst_df_2x = [pd.read_csv(f) for f in lst_files if f.startswith('ACTB  AFTER 2X')]
df_1x = pd.concat(lst_df_1x)
df_2x = pd.concat(lst_df_2x)
disp1x = df_1x.displacement_um
disp2x = df_2x.displacement_um

os.chdir(folderpath_save)
plt.figure(figsize=(9, 4), dpi=300)
sns.histplot(data=disp1x, binwidth=0.2, stat='count', color=sns.color_palette()[0], label='Isosmotic',alpha=0.5)
sns.histplot(data=disp2x, binwidth=0.2, stat='count', color=sns.color_palette()[3], label='Hyperosmotic',alpha=0.5)
plt.legend()
plt.yscale('log')
plt.title('beta-Actin CDS mRNA', fontweight='bold')
plt.xlabel('Total Displacment of an RNA Trajectory (um)')
plt.tight_layout()
plt.savefig('Total Displacement distribution.png', format='png')
