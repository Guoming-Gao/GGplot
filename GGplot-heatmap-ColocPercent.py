from os.path import dirname, join
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(color_codes=True, style="white")

path = "/Users/GGM/Documents/GGscripts/RNA-Condensate-SPT/Coloc-MaxStack/Control-Demostration/SimulatedData_ONI/Coloc_percent_200nm.csv"
coloc_threshold = 150  # nm
data = pd.read_csv(path)
data = pd.pivot(data, index="Nleft", columns="Nright")

plt.figure(dpi=300)
ticks = data.index.to_numpy(dtype=int)
ax = sns.heatmap(data=data, xticklabels=ticks, yticklabels=ticks, annot=True, fmt=".1f")
ax.invert_yaxis()
plt.xlabel("N$_{right}$ per ONI FOV", weight="bold")
plt.ylabel("N$_{left}$ per ONI FOV", weight="bold")
plt.title("<200 nm Colocalization Percentage(%)", weight="bold")
# plt.axis("scaled")
plt.tight_layout()
path_save = join(dirname(path), "Coloc_percent_heatmap_200nm.png")
plt.savefig(path_save, format="png", bbox_inches="tight")
plt.close()
