from tkinter import filedialog as fd
from os.path import join, dirname, basename
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit

sns.set(color_codes=True, style="white")

print("Type in title:")
title = input()

print("Type in the time between frames (seconds):")
t_between_frames = float(input())

print("Choose the D files to plot:")
lst_files = list(fd.askopenfilenames())

# print("Choose folder path to save:")
# folderpath = fd.askdirectory()
# os.chdir(folderpath)
folderpath = "/Volumes/AnalysisGG/PROCESSED_DATA/RNA-diffusion-in-FUS/bioFUStether-10FUS-1Mg-10Dex-RT/NoTotalRNA"
os.chdir(folderpath)

# calculate error bounds
static_err = 0.016
um_per_pxl = 0.117
link_max = 3
log10D_low = np.log10(static_err ** 2 / (4 * (t_between_frames)))
log10D_high = np.log10((um_per_pxl * link_max) ** 2 / (4 * (t_between_frames)))

# Prepare DataFrame, filter by fitting R2
lst_log10D = []
for file in lst_files:
    df_in = pd.read_csv(file)
    # R^2 filtering
    df_R2above = df_in[df_in["R2"] >= 0.7]
    lst_log10D.extend(list(df_R2above["log10D"]))

# bootstrap
dataset_size = len(lst_log10D)
lst_replicates = [
    np.random.choice(lst_log10D, round(len(lst_log10D) * 0.5), replace=False)
    for i in range(3)
]

lst_density = []
for data in lst_replicates:
    density, bins = np.histogram(
        data, bins=40, range=(log10D_low - 1.5, log10D_high + 1.5), density=True
    )
    lst_density.append(density)

array_density = np.stack(lst_density)
array_density.shape
histo_mean = np.mean(array_density, axis=0)
histo_STD = np.std(array_density, axis=0)

plt.figure(figsize=(9, 4), dpi=200)
plt.errorbar(
    x=(bins[1:] + bins[:-1]) / 2,
    y=histo_mean,
    yerr=histo_STD,
    fmt=".",
    capsize=5,
    color="dimgray",
)
plt.axvspan(
    log10D_low - 1.5, log10D_low, facecolor="dimgray", alpha=0.2, edgecolor="none"
)
plt.axvspan(
    log10D_high, log10D_high + 1.5, facecolor="dimgray", alpha=0.2, edgecolor="none"
)
plt.xlim(log10D_low - 1.5, log10D_high + 1.5)
plt.title(title, fontsize=13, fontweight="bold")
plt.xlabel("log$_{10}$D ($\mu$m^2/s)", weight="bold")
plt.tight_layout()
fsave = join(folderpath, title + ".png")
plt.savefig(fsave, format="png")
