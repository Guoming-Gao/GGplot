import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(color_codes=True, style="whitegrid")

sigma = 0.016  # unit: um
um_per_pxl = 0.117
link_max = 3  # unit: pixels

os.chdir("/Users/GGM/Documents/GGscripts/RNA-Condensate-SPT")

t_between_frames = np.arange(10, 2000, 0.1)  # unit: ms
# lower bounds determiend by static localization error
log10D_low = np.log10(sigma ** 2 / (4 * (t_between_frames / 1000)))
# higher bounds determiend by max linking length
log10D_high = np.log10((um_per_pxl * link_max) ** 2 / (4 * (t_between_frames / 1000)))
plt.figure(dpi=500)
sns.lineplot(x=t_between_frames, y=log10D_high, color="firebrick", label="Higher Bound")
sns.lineplot(x=t_between_frames, y=log10D_low, color="dimgray", label="Lower Bound")
plt.title("D Detection Limits", weight="bold")
plt.ylabel("log10D, $\mu$m$^2$/s", weight="bold")
plt.xlabel("Time Between Frames, ms", weight="bold")
plt.xticks(np.arange(0, 2100, 100), rotation=45)
plt.xlim(0, 2000)
plt.tight_layout()
plt.savefig("D-limit-freq.png", format="png")


t_between_frames = np.arange(10, 200, 0.1)  # unit: ms
# lower bounds determiend by static localization error
log10D_low = np.log10(sigma ** 2 / (4 * (t_between_frames / 1000)))
# higher bounds determiend by max linking length
log10D_high = np.log10((um_per_pxl * link_max) ** 2 / (4 * (t_between_frames / 1000)))
plt.figure(dpi=500)
sns.lineplot(x=t_between_frames, y=log10D_high, color="firebrick", label="Higher Bound")
sns.lineplot(x=t_between_frames, y=log10D_low, color="dimgray", label="Lower Bound")
plt.title("D Detection Limits", weight="bold")
plt.ylabel("log10D, $\mu$m$^2$/s", weight="bold")
plt.xlabel("Time Between Frames, ms", weight="bold")
plt.xticks(np.arange(10, 210, 10), rotation=45)
plt.xlim(0, 200)
plt.tight_layout()
plt.savefig("D-limit-freq-zoomed.png", format="png")
