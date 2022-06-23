import os
from os.path import isdir, join, dirname
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

# The code aims to compare two sets of data, it first pooling all data from each folder, then plot them together

######################################
# path and label
folder1 = "/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Analysis-20220507-PairCorr/THOR_2x"
label1 = "THOR, 300 mM Na+"
folder2 = "/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Analysis-20220507-PairCorr/L941_2x"
label2 = "L941, 300 mM Na+"

fname_save = "Cross PairCorr-Compare-THORvsL941-2x.png"

r_max_nm = 5000
dr_nm = 200
array_bin = np.arange(0, r_max_nm, dr_nm)

######################################
# pooling datasets
lst_files = [f for f in os.listdir(folder1) if f.endswith("PairCorr.csv")]
lst_array = []
for file in lst_files:
    df = pd.read_csv(join(folder1, file)).drop(columns=["roi"])
    df[df == 0] = np.nan
    lst_array.append(df.to_numpy())
array_alldata = np.concatenate(tuple(lst_array), axis=0)

array_mean1 = np.nanmean(array_alldata, axis=0)
sem1 = st.sem(df.to_numpy(), axis=0, nan_policy="omit")

lst_files = [f for f in os.listdir(folder2) if f.endswith("PairCorr.csv")]
lst_array = []
for file in lst_files:
    df = pd.read_csv(join(folder2, file)).drop(columns=["roi"])
    df[df == 0] = np.nan
    lst_array.append(df.to_numpy())
array_alldata = np.concatenate(tuple(lst_array), axis=0)

array_mean2 = np.nanmean(array_alldata, axis=0)
sem2 = st.sem(df.to_numpy(), axis=0, nan_policy="omit")

######################################
# plot
plt.figure(dpi=500)
plt.plot(array_bin, array_mean1, color="steelblue", alpha=0.9, label=label1)
plt.fill_between(
    array_bin,
    array_mean1 - sem1.data,
    array_mean1 + sem1.data,
    color="steelblue",
    alpha=0.2,
    edgecolor="none",
)
plt.plot(array_bin, array_mean2, color="firebrick", alpha=0.9, label=label2)
plt.fill_between(
    array_bin,
    array_mean2 - sem2.data,
    array_mean2 + sem2.data,
    color="firebrick",
    alpha=0.2,
    edgecolor="none",
)
plt.xlabel("r, nm")
plt.ylabel("Cross Pair Correlation, C(r)")
plt.legend()
plt.tight_layout()
path_save = join(dirname(folder1), fname_save)
plt.savefig(path_save, format="png")
plt.close()
