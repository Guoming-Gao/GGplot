import pandas as pd
import numpy as np
import os
from os.path import dirname, basename, join
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(color_codes=True, style="white")


fpath = "/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/Constructs-Primers/RNA transcripts-in vitro/GelFigs-RNAs/imageJ Gel Quant - 220815.csv"
folder_save = dirname(fpath)
os.chdir(folder_save)

data = pd.read_csv(fpath).set_index("Distance_cm")
df_plot = data[["AS3-1", "AS3-2", "AS3-3", "AS3-4"]]
fpath_save = basename(fpath).strip(".csv") + "-AS3.png"

plt.figure(figsize=(9, 4), dpi=200)
sns.lineplot(data=df_plot, dashes=False, alpha=0.7)
plt.xlabel("Distance (cm)", weight="bold")
plt.xlim(0, 6.24)
plt.legend(prop={"size": 11, "weight": "bold"})
plt.ylabel("Intensity", weight="bold")
plt.tight_layout()
plt.savefig(fpath_save, format="png")

df_plot = data[["AS1-4", "AS2-4", "AS3-4"]]
df_plot["AS1-4"] = df_plot["AS1-4"] / df_plot["AS1-4"].max()
df_plot["AS2-4"] = df_plot["AS2-4"] / df_plot["AS2-4"].max()
df_plot["AS3-4"] = df_plot["AS3-4"] / df_plot["AS3-4"].max()
fpath_save = basename(fpath).strip(".csv") + "-purity.png"

plt.figure(figsize=(9, 4), dpi=200)
sns.lineplot(data=df_plot, dashes=False, alpha=0.7)
plt.xlabel("Distance (cm)", weight="bold")
plt.xlim(0, 6.24)
plt.legend(prop={"size": 11, "weight": "bold"})
plt.ylabel("Normalized Intensity", weight="bold")
plt.tight_layout()
plt.savefig(fpath_save, format="png")
