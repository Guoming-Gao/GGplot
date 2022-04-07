import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns
sns.set_theme(style="whitegrid", palette="pastel")

path = '/Volumes/AnalysisGG/20220323_HOPSpurificationtest/forboxplot.csv'
os.chdir('/Volumes/AnalysisGG/20220323_HOPSpurificationtest/')
df = pd.read_csv(path)

plt.figure(figsize=(3, 6),dpi=200)
sns.boxplot(x='Osmolarity', y='Width', data=df, order=['150 mM Na+','450 mM Na+'])
plt.title('Condensate Diameter')
plt.ylabel('nm')
plt.tight_layout()
plt.savefig('size.png', format='png')
plt.close()

plt.figure(figsize=(3, 6),dpi=200)
sns.boxplot(x='Osmolarity', y='Intens', data=df, order=['150 mM Na+','450 mM Na+'])
plt.title('Condensate Intensity')
plt.ylabel('A.U.')
plt.tight_layout()
plt.savefig('intensity.png', format='png')
plt.close()
