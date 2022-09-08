import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns

sns.set(color_codes=True, style="white")

###############################################
# Loading Data
R2threshold = 0.9

folderpath = "/Volumes/AnalysisGG/PROCESSED_DATA/DiffusionCoeffcientLimit"
os.chdir(folderpath)
labels = ["2 s", "100 ms", "40 ms", "20 ms"]
fname_starts = ["FL_2s", "FL_100ms", "FL_40ms", "FL_20ms"]

lst_dataclass = []
lst_log10D = []
for label, start in zip(labels, fname_starts):
    lst_df = [
        pd.read_csv(f)
        for f in os.listdir(folderpath)
        if f.endswith("linregress_D.csv") & f.startswith(start)
    ]
    for df in lst_df:
        df_R2above = df[df["R2"] >= R2threshold]
        lst_log10D.extend(list(df_R2above["log10D"]))
        lst_dataclass.extend(list(np.repeat(label, df_R2above.shape[0])))

df_plot = pd.DataFrame({"log10D": lst_log10D, "dataclass": lst_dataclass},)
# Count number of molecules and put in label
final_labels = []
for label in labels:
    N = df_plot[df_plot["dataclass"] == label].shape[0]
    final_labels.extend(list(np.repeat(label + ", N=" + str(N), N)))
df_plot["Exposure Time"] = final_labels
df_plot = df_plot.astype(dtype={"log10D": "float64", "Exposure Time": "string"})

plt.figure(figsize=(9, 4), dpi=200)
sns.histplot(
    data=df_plot,
    x="log10D",
    hue="Exposure Time",
    common_norm=False,
    stat="probability",
)
# Theoretical Thresholds
lst_exptime = [2, 0.1, 0.04, 0.02]
for exptime, idx in zip(lst_exptime, range(4)):
    bound = np.log10((0.016 ** 2) / (4 * exptime))
    plt.axvline(x=bound, color=sns.color_palette()[idx], ls="--")
plt.title("D Lower Bounds by Imaging Frequency", weight="bold")
plt.xlabel("log10D ($\mu$m^2/s)", weight="bold")
plt.tight_layout()
plt.savefig("D lower bounds by imaging frequency.png", format="png")
