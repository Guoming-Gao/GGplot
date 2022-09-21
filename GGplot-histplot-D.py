from tkinter import filedialog as fd
from os.path import join, dirname
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set(color_codes=True, style="white")

print("Type in the title:")
title = input()

print("Do you want to set range by detection limits? (Y/N)")
set_range = input() == "Y"

if set_range:
    print("Type in the time between frames (seconds):")
    t_between_frames = float(input())

print("Choose the D files to plot:")
lst_files = list(fd.askopenfilenames())

# Prepare DataFrame, filter by fitting R2
lst_log10D = []
for file in lst_files:
    df_in = pd.read_csv(file)
    # R^2 filtering
    df_R2above = df_in[df_in["R2"] >= 0.7]
    lst_log10D.extend(list(df_R2above["log10D"]))

df_plot = pd.DataFrame({"log10D (um^2/s)": lst_log10D,}, dtype=float)


######################################
# Total D distribution with fitting
plt.figure(figsize=(9, 4), dpi=200)
if set_range:
    sigma = 0.016
    um_per_pxl = 0.117
    link_max = 3
    log10D_low = np.log10(sigma ** 2 / (4 * (t_between_frames)))
    log10D_high = np.log10((um_per_pxl * link_max) ** 2 / (4 * (t_between_frames)))
    g = sns.histplot(
        data=df_plot,
        x="log10D (um^2/s)",
        fill=True,
        stat="count",
        alpha=0.7,
        color="dimgray",
        bins=30,
        binrange=(log10D_low, log10D_high),
    )
    plt.xlim(log10D_low, log10D_high)
else:
    g = sns.histplot(
        data=df_plot,
        x="log10D (um^2/s)",
        fill=True,
        stat="count",
        alpha=0.7,
        color="dimgray",
        bins=30,
    )
plt.title(title, fontsize=13, fontweight="bold")
plt.xlabel("log10D ($\mu$m^2/s)", weight="bold")
plt.tight_layout()
fsave = join(dirname(lst_files[0]), title + ".png")
plt.savefig(fsave, format="png")
