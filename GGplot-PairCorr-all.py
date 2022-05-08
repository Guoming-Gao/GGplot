import os
from os.path import isdir, join
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from rich.progress import track

# The code aims to plot all pair correlation in all subfolders (each represent a dataset to be pooled)

path = "/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/Analysis-20220507-PairCorr/"

lst_folders = [f for f in os.listdir(path) if isdir(join(path, f))]
array_bin = np.arange(0, 5000, 200)

folder = lst_folders[0]
for folder in lst_folders:
    lst_files = [
        f for f in os.listdir(join(path, folder)) if f.endswith("PairCorr.csv")
    ]

    lst_array = []
    for file in track(lst_files, description=folder):
        df = pd.read_csv(join(path, folder, file)).drop(columns=["roi"])
        df[df == 0] = np.nan
        lst_array.append(df.to_numpy())

    array_alldata = np.concatenate(tuple(lst_array), axis=0)

    array_mean = np.nanmean(array_alldata, axis=0)
    sem = st.sem(df.to_numpy(), axis=0, nan_policy="omit")

    # saving plot data
    df_plotdata = pd.DataFrame(columns=["bins", "mean", "SEM"])
    df_plotdata["bins"] = array_bin
    df_plotdata["mean"] = array_mean
    df_plotdata["SEM"] = sem.data
    path_save = join(path, folder, folder + "-Cross-PairCorr-toplot.csv")
    df_plotdata.to_csv(path_save, index=False)

    plt.figure(dpi=500)
    plt.plot(array_bin, array_mean, color="firebrick")
    plt.fill_between(
        array_bin, array_mean - sem.data, array_mean + sem.data, color="gray", alpha=0.2
    )
    plt.title(folder)
    plt.xlabel("r, nm")
    plt.ylabel("g(r)")
    plt.tight_layout()
    path_save = join(path, folder, folder + "-Cross-PairCorr.png")
    plt.savefig(path_save, format="png")
    plt.close()
