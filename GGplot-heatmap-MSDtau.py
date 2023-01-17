import os
import matplotlib.pyplot as plt
from tkinter import filedialog as fd
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.colors as colors

sns.set(color_codes=True, style="white")

print("Select a MSD-tau-alltracks.csv file to be plot")
fpath = fd.askopenfilename()
# fpath = "/Volumes/nwalter-group/Guoming Gao/PROCESSED_DATA/RNA-diffusion-in-FUS/bioFUStether-10FUS-1Mg-10Dex-RT/NoTotalRNA/FL-0min/tracks-rawRNA-minlength-5/MSD-tau-alltracks.csv"
print("Type in a title below:")
title = input()

data = pd.read_csv(fpath)

histogram_by_tau = []
bin_centers_by_tau = []
for tau in np.arange(1, data.tau.max() + 1):
    all_MSD_per_tau = np.array(data[data.tau == tau].MSD_um2)
    hist, bin_edges = np.histogram(
        all_MSD_per_tau,
        bins=int(data.tau.max()),
        range=(data.MSD_um2.min(), data.MSD_um2.max()),
    )
    histogram_by_tau.append(hist)
    bin_centers_by_tau.append((bin_edges[1:] + bin_edges[:-1]) / 2)
matrix_MSD = np.flipud(np.stack(histogram_by_tau, axis=1))
matrix_MSD_centers = np.flipud(np.stack(bin_centers_by_tau, axis=1))
fig, ax = plt.subplots(dpi=600)
xmin = 0
xmax = data.tau.max()
ymin = data.MSD_um2.min()
ymax = data.MSD_um2.max()
pc = plt.imshow(
    matrix_MSD,
    cmap=plt.cm.Reds,
    extent=[xmin, xmax, ymin, ymax],
    aspect="auto",
    norm=colors.LogNorm(vmin=1, vmax=matrix_MSD.max() / 5e2),
)
fig.colorbar(pc, ax=ax)
plt.title(title, weight="bold")
plt.xlabel("tau, frames", weight="bold")
plt.ylabel("MSD, $\mu$m^2", weight="bold")
plt.tight_layout()
fpath_save = fpath.strip(".csv") + "-heatmap.png"
plt.savefig(fpath_save, format="png")
