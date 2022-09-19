import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

sns.set(color_codes=True, style="white")

title = "RNA in FUS\n-low frequency"
folderpath = "/Volumes/AnalysisGG/PROCESSED_DATA/Static-RNA-Control/StaticPercentEstimate/example-FUSinCondensate/"
os.chdir(folderpath)
lst_files = [
    f
    for f in os.listdir(folderpath)
    if f.endswith("linregress_D.csv") & f.startswith("low")
]

lst_tot = []
lst_immobile = []
for f in lst_files:
    df = pd.read_csv(f)
    lst_tot.append(df.shape[0])
    df_slope = df[df["slope"] > 0]
    df_R2 = df_slope[df_slope["R2"] > 0.7]
    lst_immobile.append(df.shape[0] - df_R2.shape[0])

ratio = 100 * np.array(lst_immobile) / np.array(lst_tot)
data = pd.DataFrame({"percent": ratio}, dtype=float)
plt.figure(figsize=(3, 6), dpi=200)
sns.barplot(data=data, y="percent", alpha=0.7, color="dimgray", capsize=0.4)
sns.stripplot(data=data, y="percent", color="dimgray")
plt.ylabel("Immobile Molecules %", weight="bold")
plt.ylim(0, 100)
plt.title(title, weight="bold")
plt.xlabel("")
plt.tight_layout()
plt.savefig("ImmobilePercent.png", format="png")
