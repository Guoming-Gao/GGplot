from tkinter import filedialog as fd
from os.path import join, dirname, basename
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set(color_codes=True, style="white")

# lst_files = list(fd.askopenfilenames())
file = "/Volumes/AnalysisGG/PROCESSED_DATA/CondensateAnalysis/Dcp1a-2x-2s/20221015-UGD-100msexposure-2sperframe-FOV-1.csv"

df = pd.read_csv(file)
lst_meanInt = df.amp
lst_estD = df.sig
lst_t = df.t
# for file in lst_files:
# timetag = basename(file).strip(".csv").split("_")[-1]
# df = pd.read_csv(file)
# lst_meanInt.extend(list(df["meanIntensity"]))
# lst_estD.extend(list(df["estDiameter"]))
# lst_t.extend(list((df["t"] + int(timetag) * 200) * 2 / 60))  # unit: min

data = pd.DataFrame({"t": lst_t, "meanInt": lst_meanInt, "estD": lst_estD})

lst_estD_prob = []
lst_meanInt_prob = []
timepoints = np.sort(data.t.unique())
for t in timepoints:
    current_frame = data[data.t == t]
    estD_density, estD_bins = np.histogram(
        current_frame["estD"], bins=10, range=(2, 5), density=True
    )
    lst_estD_prob.append(estD_density * (5 - 2) / 10)
    meanInt_density, meanInt_bins = np.histogram(
        current_frame["meanInt"], bins=10, range=(400, 1200), density=True
    )
    lst_meanInt_prob.append(meanInt_density * (1200 - 400) / 10)


array_meanInt_prob = np.stack(lst_meanInt_prob)
df_array = pd.DataFrame(
    np.transpose(array_meanInt_prob),
    columns=[round(x) for x in timepoints],
    index=[round(x) for x in (meanInt_bins[1:] + meanInt_bins[:-1]) / 2],
)
plt.figure(1)
ax = sns.heatmap(df_array, robust=True)
ax.invert_yaxis()
plt.ylabel("Condensate Intensity", weight="bold")
plt.xlabel("Time (min)", weight="bold")
plt.savefig(join(dirname(file), "heatmap-condensate intensity.png"), format="png")

array_estD_prob = np.stack(lst_estD_prob)
df_array = pd.DataFrame(
    np.transpose(array_estD_prob),
    columns=[round(x) for x in timepoints],
    index=[round(x) for x in (estD_bins[1:] + estD_bins[:-1]) / 2 * 117],
)
plt.figure(2)
ax = sns.heatmap(df_array, robust=True)
ax.invert_yaxis()
plt.ylabel("Estimated Diameter (nm)", weight="bold")
plt.xlabel("Time (min)", weight="bold")
plt.savefig(join(dirname(file), "heatmap-condensate estD.png"), format="png")

lst_N = [data[data.t == t].shape[0] for t in timepoints]
lst_condensate_Int_median = [data[data.t == t]["meanInt"].mean() for t in timepoints]
lst_condensate_estD_median = [
    data[data.t == t]["estD"].mean() * 0.117 for t in timepoints
]
plt.figure(3)
plt.plot(timepoints, lst_N)
plt.xlabel("Time (min)", weight="bold")
plt.ylabel("Number of Condensates", weight="bold")
plt.savefig(join(dirname(file), "line-condensate N.png"), format="png")

plt.figure(4)
plt.plot(timepoints, lst_condensate_Int_median)
plt.xlabel("Time (min)", weight="bold")
plt.ylabel("Median Condensate Intensity (pixel)", weight="bold")
plt.savefig(join(dirname(file), "line-condensate intensity.png"), format="png")

plt.figure(5)
plt.plot(timepoints, lst_condensate_estD_median)
plt.xlabel("Time (min)", weight="bold")
plt.ylabel("Median Estimated Diameter (pixel)", weight="bold")
plt.savefig(join(dirname(file), "line-condensate estD.png"), format="png")

plt.figure(6)
plt.plot(
    timepoints,
    np.array(lst_N) * np.pi * (np.array(lst_condensate_estD_median) / 2) ** 2,
)
plt.xlabel("Time (min)", weight="bold")
plt.ylabel("Median Condensate Area ($\mu$m$^2$)", weight="bold")
plt.savefig(join(dirname(file), "line-condensate area.png"), format="png")
