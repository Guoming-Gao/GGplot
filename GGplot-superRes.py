import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from progress.bar import IncrementalBar

folderpath_save = '/Volumes/nwalter-group/Shankar/13. Exosome project/Exosome data/Guoming_processing/negative_results'
os.chdir(folderpath_save)

fpath = '/Volumes/nwalter-group/Shankar/13. Exosome project/Exosome data/Guoming_processing/negative.csv'
figsize = np.array([428,684], dtype=float)

df = pd.read_csv(fpath)
x = df.POSITION_X.to_numpy(dtype=float)
y = df.POSITION_Y.to_numpy(dtype=float)

# plot the general overview
plt.figure(num=1, figsize=(4,7), dpi=300)
plt.scatter(x, y, c='firebrick', alpha=0.5, s=0.01, edgecolors='none')
plt.axis('scaled')
plt.xlim(0, figsize[0])
plt.ylim(0, figsize[1])
plt.savefig('general_overview.png', format='png')

# plot detail images
plt.figure(num=2, figsize=(5,5), dpi=300)
plt.scatter(x, y, c='firebrick', alpha=0.5, s=1, edgecolors='none')
plt.axis('scaled')
range_x = np.arange(0, figsize[0], 20, dtype=int)
range_y = np.arange(0, figsize[1], 20, dtype=int)
max = (len(range_x)-1) * (len(range_y)-1)
bar = IncrementalBar('Plotting detailed images...', suffix='%(percent).1f%% - %(eta)ds', max=max)
N = 0
for i in range(len(range_x)-1):
    for j in range(len(range_y)-1):
        plt.xlim(range_x[i], range_x[i+1])
        plt.ylim(range_y[j], range_y[j+1])
        fsave = str(N) + '.png'
        plt.savefig(fsave, format='png')
        N += 1
        bar.next()
bar.finish()
