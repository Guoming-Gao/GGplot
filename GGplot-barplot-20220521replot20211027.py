import os
from os.path import join, basename, dirname
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import seaborn as sns

sns.set(color_codes=True, style="white")


path = "/Volumes/AnalysisGG/PROCESSED DATA/2021Nov-LiveCell-WhyIsWrong/ALL-Coloc-Results.csv"

df_all = pd.read_csv(path)
# Omit Xin's L941 data:
df_all = df_all[df_all.RNAtype != "L941 lncRNA"]


data = df_all[df_all.HOPScondition == "300 mM Na+"]

# dict = {
# "ActB": "beta-Actin mRNA",
# "L941": "L941 lncRNA",
# "THOR-d": "THOR-delta lncRNA",
# "THOR": "THOR lncRNA",
# "miR-21gs": "miR-21 guide strand",
# "miR-21ds": "miR-21 double strand",
# "SOX2": "SOX2 mRNA",
# }
# data.loc[:, "RNAtype"] = data.RNAtype.replace(dict)

os.chdir(dirname(path))
plt.figure(figsize=(5, 6), dpi=200)
sns.barplot(
    data=data,
    x="RNAtype",
    y="N_total in roi",
    order=["miR-21ds", "miR-21gs", "THOR", "THOR-d", "L941", "ActB", "SOX2"],
)
sns.stripplot(
    data=data,
    x="RNAtype",
    y="N_total in roi",
    order=["miR-21ds", "miR-21gs", "THOR", "THOR-d", "L941", "ActB", "SOX2"],
    color=".25",
)
# plt.yscale("log")
plt.ylim(0, 500)
plt.xticks(rotation=30)
plt.title("Number of RNA per cell")
plt.tight_layout()

plt.savefig("N RNA per cell-foranalysis20211027.svg", format="svg")
