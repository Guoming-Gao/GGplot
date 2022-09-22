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

print("Type in the time between frames (seconds):")
t_between_frames = float(input())

print("Choose the D files to plot:")
lst_files = list(fd.askopenfilenames())


def Gauss(x, A, x0, sigma):
    return A * np.exp((-1 / 2) * ((x - x0) ** 2 / sigma ** 2))


def DualGauss(x, A1, x1, sigma1, A2, x2, sigma2):
    return A1 * np.exp((-1 / 2) * ((x - x1) ** 2 / sigma1 ** 2)) + A2 * np.exp(
        (-1 / 2) * ((x - x2) ** 2 / sigma2 ** 2)
    )


def DualGauss_fit_plot_text(data, range):
    counts, bins = np.histogram(data, bins=40, range=range)
    (A1, x1, sigma1, A2, x2, sigma2), pcov = curve_fit(
        DualGauss, (bins[1:] + bins[:-1]) / 2, counts
    )
    err_A1, err_x1, err_sigma1, err_A2, err_x2, err_sigma2 = np.sqrt(np.diag(pcov))
    curve_x = np.arange(bins[0], bins[-1], 0.01)
    curve_ydual = DualGauss(curve_x, A1, x1, sigma1, A2, x2, sigma2)
    curve_y1 = Gauss(curve_x, A1, x1, sigma1)
    curve_y2 = Gauss(curve_x, A2, x2, sigma2)
    plt.plot(curve_x, curve_ydual, color="dimgray", linewidth=2)
    plt.plot(curve_x, curve_y1, color=sns.color_palette()[1], linewidth=2)
    plt.plot(curve_x, curve_y2, color=sns.color_palette()[3], linewidth=2)
    plt.axvline(x=x1, color=sns.color_palette()[1], ls="--")
    plt.axvline(x=x2, color=sns.color_palette()[3], ls="--")

    # label with text
    plt.text(
        0.76,
        0.83,
        "log10D$_1$ = " + str(round(x1, 2)) + "$\pm$" + str(round(err_x1, 2)),
        weight="bold",
        fontsize=15,
        color=sns.color_palette()[1],
        transform=plt.gcf().transFigure,
    )
    plt.text(
        0.76,
        0.75,
        "log10D$_2$ = " + str(round(x2, 2)) + "$\pm$" + str(round(err_x2, 2)),
        weight="bold",
        fontsize=15,
        color=sns.color_palette()[3],
        transform=plt.gcf().transFigure,
    )


def Gauss_fit_plot_text(data, range):
    counts, bins = np.histogram(data, bins=40, range=range)
    # Fit to Gauss and plot individually and combined
    (A, x0, sigma), pcov = curve_fit(Gauss, (bins[1:] + bins[:-1]) / 2, counts)
    err_A, err_x0, err_sigma = np.sqrt(np.diag(pcov))
    curve_x = np.arange(bins[0], bins[-1], 0.01)
    curve_y = Gauss(curve_x, A, x0, sigma)
    plt.plot(curve_x, curve_y, color=sns.color_palette()[3], linewidth=2)
    plt.axvline(x=x0, color=sns.color_palette()[3], ls="--")

    # label with text
    plt.text(
        0.76,
        0.83,
        "$\mu$ = " + str(round(x0, 2)) + "$\pm$" + str(round(err_x0, 2)),
        weight="bold",
        fontsize=15,
        color=sns.color_palette()[3],
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


# calculate error bounds
static_err = 0.016
um_per_pxl = 0.117
link_max = 3
log10D_low = np.log10(static_err ** 2 / (4 * (t_between_frames)))
log10D_high = np.log10((um_per_pxl * link_max) ** 2 / (4 * (t_between_frames)))

# estimate if dual Gaussian fitting is better than single Gaussian
counts, bins = np.histogram(
    df_plot["log10D (um^2/s)"], bins=40, range=(log10D_low - 1.5, log10D_high + 1.5),
)
(A1, x1, sigma1, A2, x2, sigma2), pcov = curve_fit(
    DualGauss, (bins[1:] + bins[:-1]) / 2, counts
)
if abs(x1 - x2) > 0.5:
    DualFit = True
else:
    DualFit = False

plt.figure(figsize=(9, 4), dpi=200)
g = sns.histplot(
    data=df_plot,
    x="log10D (um^2/s)",
    fill=True,
    stat="count",
    alpha=0.5,
    color="dimgray",
    bins=40,
    binrange=(log10D_low - 1.5, log10D_high + 1.5),
)
plt.axvspan(log10D_low - 1.5, log10D_low, color="dimgray", alpha=0.3)
plt.axvspan(log10D_high, log10D_high + 1.5, color="dimgray", alpha=0.3)
plt.xlim(log10D_low - 1.5, log10D_high + 1.5)
if DualFit:
    DualGauss_fit_plot_text(
        data=df_plot["log10D (um^2/s)"], range=(log10D_low - 1.5, log10D_high + 1.5)
    )
else:
    Gauss_fit_plot_text(
        data=df_plot["log10D (um^2/s)"], range=(log10D_low - 1.5, log10D_high + 1.5),
    )
plt.title(title, fontsize=13, fontweight="bold")
plt.xlabel("log10D ($\mu$m^2/s)", weight="bold")
plt.tight_layout()
fsave = join(dirname(lst_files[0]), title + ".png")
plt.savefig(fsave, format="png")
