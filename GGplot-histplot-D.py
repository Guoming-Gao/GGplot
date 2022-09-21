from tkinter import filedialog as fd
from os.path import join, dirname
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit

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


def Gauss_fit_plot_text(data, range, color, text_x, text_y):
    try:
        counts, bins = np.histogram(data, bins=30, range=range)
    except:
        counts, bins = np.histogram(data, bins=30)

    def Gauss(x, A, x0, sigma):
        return A * np.exp((-1 / 2) * ((x - x0) ** 2 / sigma ** 2))

    # Fit to Gauss and plot individually and combined
    (A, x0, sigma), pcov = curve_fit(Gauss, (bins[1:] + bins[:-1]) / 2, counts,)
    err_A, err_x0, err_sigma = np.sqrt(np.diag(pcov))
    curve_x = np.arange(bins[0], bins[-1], 0.01)
    curve_y = Gauss(curve_x, A, x0, sigma)
    plt.plot(curve_x, curve_y, color=color, linewidth=2)
    plt.axvline(x=x0, color=color, ls="--")

    # label with text
    plt.text(
        text_x,
        text_y,
        "$\mu$ = " + str(round(x0, 2)) + "$\pm$" + str(round(err_x0, 2)),
        weight="bold",
        fontsize=11,
        color=color,
        transform=plt.gcf().transFigure,
    )


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
    static_err = 0.016
    um_per_pxl = 0.117
    link_max = 3
    log10D_low = np.log10(static_err ** 2 / (4 * (t_between_frames)))
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
    Gauss_fit_plot_text(
        data=df_plot["log10D (um^2/s)"],
        range=(log10D_low, log10D_high),
        color="firebrick",
        text_x=0.815,
        text_y=0.825,
    )
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
    data = df_plot["log10D (um^2/s)"]
    Gauss_fit_plot_text(data, None, "firebrick", 0.815, 0.825)
plt.title(title, fontsize=13, fontweight="bold")
plt.xlabel("log10D ($\mu$m^2/s)", weight="bold")
plt.tight_layout()
fsave = join(dirname(lst_files[0]), title + ".png")
plt.savefig(fsave, format="png")
