import os
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from rich.progress import track
from scipy.optimize import curve_fit

sns.set(color_codes=True, style="white")

folderpath = "/Volumes/AnalysisGG/PROCESSED_DATA/RNA-diffusion-in-FUS/20220712_FLmRNA_10FUS_1Mg_10Dex_noTotR_24C/"
os.chdir(folderpath)
lst_fname = [f for f in os.listdir(folderpath) if f.endswith("RNA_linregress_D.csv")]
R2threshold = 0.9

######################################
# Prepare DataFrame, filter by fitting R2
lst_dataclass = []
lst_log10D = []
for fname in lst_fname:
    if fname.startswith("high"):
        dataclass = "50 Hz"
    elif fname.startswith("low"):
        dataclass = "0.5 Hz"
    df_in = pd.read_csv(join(folderpath, fname))
    # R^2 filtering
    df_R2above = df_in[df_in["R2"] >= R2threshold]
    lst_dataclass.extend(list(np.repeat(dataclass, df_R2above.shape[0])))
    lst_log10D.extend(list(df_R2above["log10D"]))

df_out = pd.DataFrame({"dataclass": lst_dataclass, "log10D (um^2/s)": lst_log10D,},)

dfslow = df_out[df_out["dataclass"] == "0.5 Hz"]
dffast = df_out[df_out["dataclass"] == "50 Hz"]
Nslow = dfslow.shape[0]
Nfast = dffast.shape[0]
tag_slow = "0.5 Hz, N=" + str(Nslow)
tag_fast = "50 Hz, N=" + str(Nfast)
df_plot = pd.DataFrame(
    {
        "Frequency": list(np.repeat(tag_slow, Nslow))
        + list(np.repeat(tag_fast, Nfast)),
        "log10D (um^2/s)": list(dfslow["log10D (um^2/s)"])
        + list(dffast["log10D (um^2/s)"]),
    },
)


######################################
# Total D distribution with fitting
plt.figure(figsize=(9, 4), dpi=200)
g = sns.histplot(
    data=df_plot,
    x="log10D (um^2/s)",
    hue="Frequency",
    fill=True,
    stat="Proba",
    alpha=0.7,
    common_norm=True,
    bins=30,
)
plt.title(
    "Diffusion Coefficient Distribution \n(only within condensates)",
    fontsize=13,
    fontweight="bold",
)

# calculate limits
# lower bounds determiend by static localization error 55 nm
low1 = np.log10(0.055 ** 2 / (4 * 2))
low2 = np.log10(0.055 ** 2 / (4 * 0.02))
# higher bounds determiend by max linking length 3 pixels
high1 = np.log10((0.117 * 3) ** 2 / (4 * 2))
high2 = np.log10((0.117 * 3) ** 2 / (4 * 0.02))
# plot limits
plt.axvline(x=low1, color="dimgray", ls=":")
plt.axvline(x=low2, color="dimgray", ls=":")
plt.axvline(x=high1, color="firebrick", ls=":")
plt.axvline(x=high2, color="firebrick", ls=":")

plt.xlabel("log10D ($\mu$m^2/s)", weight="bold")
# g.legend_.set_title(None)
plt.tight_layout()
fsave = "D Distribution-dual frequency.png"
plt.savefig(fsave, format="png")
