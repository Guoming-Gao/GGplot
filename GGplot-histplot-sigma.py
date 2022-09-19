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

lst_sigma = []
for f in lst_files:
    df = pd.read_csv(f)
    df_slope = df[df["slope"] > 0]
    df_R2 = df_slope[df_slope["R2"] > 0.7]
    lst_sigma.extend(list(df_R2["sigma(nm)"]))

plt.figure(figsize=(9, 4), dpi=200)
plt.hist(
    lst_sigma, bins=30, color="dimgray", alpha=0.7, density=False,
)
plt.axvline(x=16, color="firebrick", ls=":", alpha=0.7)
plt.title(r"MSD-$\tau$ Fitting Derived $\sigma$, low", weight="bold")
plt.xlabel("$\sigma$ (nm)", weight="bold")
plt.ylabel("Counts", weight="bold")
plt.tight_layout()
plt.savefig("MSD-tau Fitting Derived Sigma-low.png", format="png")
