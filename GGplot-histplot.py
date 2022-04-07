import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns
sns.set(color_codes=True, style='white')

###############################################
# Loading Data
folderpath = '/Volumes/AnalysisGG/KNIME_results_from_greatlakes/20210909-Condensates'
os.chdir(folderpath)
lst_files = [f for f in os.listdir(folderpath) if f.endswith("filtered.csv")]


data = np.array([], dtype=float)
for file in lst_files:
    df = pd.read_csv(file)
    data = np.append(data, df.meanIntensity.to_numpy(dtype=float))

plt.figure(figsize=(9, 4), dpi=200)
sns.histplot(data=data, binwidth=100, stat='probability', color=sns.color_palette()[1])
plt.title('Condensates')
plt.xlim(0,3000)
plt.xlabel('Mean Intensity, A.U.')
plt.tight_layout()
plt.savefig('Condensates mean intensity histogram.png', format='png')
