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
lst_files = [f for f in os.listdir(folderpath) if f.endswith("_angles.csv")]

lst_df_1x = [pd.read_csv(f) for f in lst_files if f.startswith('ACTB  BEFORE-1X')]
lst_df_2x = [pd.read_csv(f) for f in lst_files if f.startswith('ACTB  AFTER 2X')]
df_1x = pd.concat(lst_df_1x)
df_2x = pd.concat(lst_df_2x)

# angles = df_1x[(df_1x.displacement_um<0.9)&(df_1x.displacement_um>0.8)].angles

os.chdir(folderpath_save)

angles = df_1x[df_1x.displacement_um>2].angles
hist, bin_edges = np.histogram(angles, bins=36, range=(-180,180), density=True)
plt.figure(figsize=(5, 5), dpi=300)
N = 36
theta = np.linspace(-np.pi, np.pi, N, endpoint=False)
ax = plt.subplot(projection='polar')
ax.bar(theta, hist, width=0.17, bottom=0.001, color='black', alpha=0.5)
ax.set_yticks([])
plt.ylim(0,0.0065)
plt.tight_layout()
plt.savefig('Angle distribution_1x_2up.png', format='png')
plt.close()
print('1x >2', angles.size)



angles = df_1x[df_1x.displacement_um<=2].angles
hist, bin_edges = np.histogram(angles, bins=36, range=(-180,180), density=True)
plt.figure(figsize=(5, 5), dpi=300)
N = 36
theta = np.linspace(-np.pi, np.pi, N, endpoint=False)
ax = plt.subplot(projection='polar')
ax.bar(theta, hist, width=0.17, bottom=0.001, color='black', alpha=0.5)
ax.set_yticks([])
plt.ylim(0,0.0065)
plt.tight_layout()
plt.savefig('Angle distribution_1x_2down.png', format='png')
plt.close()
print('1x <2', angles.size)


angles = df_2x[df_2x.displacement_um>2].angles
hist, bin_edges = np.histogram(angles, bins=36, range=(-180,180), density=True)
plt.figure(figsize=(5, 5), dpi=300)
N = 36
theta = np.linspace(-np.pi, np.pi, N, endpoint=False)
ax = plt.subplot(projection='polar')
ax.bar(theta, hist, width=0.17, bottom=0.001, color='black', alpha=0.5)
ax.set_yticks([])
plt.ylim(0,0.0065)
plt.tight_layout()
plt.savefig('Angle distribution_2x_2up.png', format='png')
plt.close()
print('2x >2', angles.size)


angles = df_2x[df_2x.displacement_um<=2].angles
hist, bin_edges = np.histogram(angles, bins=36, range=(-180,180), density=True)
plt.figure(figsize=(5, 5), dpi=300)
N = 36
theta = np.linspace(-np.pi, np.pi, N, endpoint=False)
ax = plt.subplot(projection='polar')
ax.bar(theta, hist, width=0.17, bottom=0.001, color='black', alpha=0.5)
ax.set_yticks([])
plt.ylim(0,0.0065)
plt.tight_layout()
plt.savefig('Angle distribution_2x_2down.png', format='png')
plt.close()
print('2x <2', angles.size)


plt.figure(figsize=(9, 5), dpi=300)
sns.histplot(df_1x.angles, bins=36, binrange=(-180,180), stat='probability', color=sns.color_palette()[0], label='Isosmotic',alpha=0.5)
sns.histplot(df_2x.angles, bins=36, binrange=(-180,180), stat='probability', color=sns.color_palette()[3], label='Hyperosmotic',alpha=0.5)
plt.legend()
plt.yscale('log')
plt.xticks(np.arange(-180,181,30))
plt.xlim(-180,180)
plt.title('beta-Actin CDS mRNA', fontweight='bold')
plt.xlabel('Angle between steps (degree)')
plt.tight_layout()
plt.savefig('Angle distribution_1xallvs2xall.png', format='png')
