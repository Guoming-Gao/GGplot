import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(color_codes=True, style="whitegrid")

folder = "/Volumes/nwalter-group/Guoming Gao/PROCESSED_DATA/RNA-diffusion-in-FUS/bioFUStether-10FUS-1Mg-10Dex-RT/NoTotalRNA/"
os.chdir(folder)
lst_fname = [
    "FL-saSPT_output_0hr.csv",
    "FL-saSPT_output_1hr.csv",
    "FL-saSPT_output_2hr.csv",
    "FL-saSPT_output_3hr.csv",
]
lst_tag = [
    "0 hr",
    "1 hr",
    "2 hr",
    "3 hr",
]
title = "saSPT-FL-10Dex-noTotalRNA-aging"

plt.figure(dpi=500)
for i in range(len(lst_fname)):
    fname = lst_fname[i]
    tag = lst_tag[i]
    df = pd.read_csv(fname)
    log10D = df["log10D"].to_numpy(dtype=float)
    prob = df["Density"].to_numpy(dtype=float)

    sns.lineplot(
        x=log10D, y=prob, color=sns.color_palette()[i], label=tag, alpha=0.5, lw=3
    )

plt.yscale("log")
plt.title(title, weight="bold")
plt.ylabel("Probability", weight="bold")
plt.xlabel("log10D, $\mu$m$^2$/s", weight="bold")
plt.tight_layout()
plt.savefig(title + ".png", format="png")
