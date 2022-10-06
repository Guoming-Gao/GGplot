import os
from os.path import dirname
from tkinter import filedialog as fd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns

sns.set(color_codes=True, style="white")

###############################################
# Loading Data

print("Type in Dataset Name:")
dataset_name = input()

print("Choose the D files to plot:")
lst_files = list(fd.askopenfilenames())

print("Choose folder path to save:")
folderpath = fd.askdirectory()
os.chdir(folderpath)

lst_slope = []
lst_R2 = []
for f in lst_files:
    df = pd.read_csv(f)
    lst_slope.extend(list(df["slope"]))
    lst_R2.extend(list(df["R2"]))

plt.figure(figsize=(9, 4), dpi=200)
plt.hist(
    lst_slope,
    bins=np.arange(-0.025, 0.025, 0.0005),
    color="dimgray",
    alpha=0.7,
    density=False,
)
plt.axvline(x=0, color="firebrick", ls=":", alpha=0.7)
plt.title(r"MSD-$\tau$ Fitting Slope", weight="bold")
plt.xlabel("slope ('$\mu$m^2/s')", weight="bold")
plt.ylabel("Counts", weight="bold")
plt.tight_layout()
plt.savefig(dataset_name + "_MSD-tau Fitting Slope.png", format="png")


plt.figure(figsize=(9, 4), dpi=200)
plt.hist(lst_R2, bins=50, color="dimgray", alpha=0.7, density=False)
plt.axvline(x=0.7, color="firebrick", ls=":", alpha=0.7)
plt.title(r"MSD-$\tau$ Fitting R$^2$", weight="bold")
plt.xlabel("R$^2$", weight="bold")
plt.ylabel("Counts", weight="bold")
plt.tight_layout()
plt.savefig(dataset_name + "_MSD-tau Fitting R2.png", format="png")
