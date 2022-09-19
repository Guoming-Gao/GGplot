import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns

sns.set(color_codes=True, style="white")

###############################################
# Loading Data
folderpath = "/Volumes/AnalysisGG/PROCESSED_DATA/Static-RNA-Control/StaticPercentEstimate/example-FUSinCondensate/"
os.chdir(folderpath)
lst_files = [
    f
    for f in os.listdir(folderpath)
    if f.endswith("linregress_D.csv") & f.startswith("low")
]

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
plt.savefig("MSD-tau Fitting Slope.png", format="png")


plt.figure(figsize=(9, 4), dpi=200)
plt.hist(lst_R2, bins=50, color="dimgray", alpha=0.7, density=False)
plt.axvline(x=0.7, color="firebrick", ls=":", alpha=0.7)
plt.title(r"MSD-$\tau$ Fitting R$^2$", weight="bold")
plt.xlabel("R$^2$", weight="bold")
plt.ylabel("Counts", weight="bold")
plt.tight_layout()
plt.savefig("MSD-tau Fitting R2.png", format="png")
