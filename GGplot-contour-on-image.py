from tkinter import filedialog as fd
from os.path import join, dirname, basename
from tifffile import imread
import matplotlib.pyplot as plt
from skimage import exposure
import numpy as np
import cv2
import pickle

path_contours = "/Volumes/AnalysisGG/PROCESSED_DATA/2022July-RNAinFUS-preliminary/20220712_FLmRNA_10FUS_1Mg_10Dex_noTotR_24C/high_freq_50Hz_both_on_FOV-5-condensate_contours.pkl"

contours, img = pickle.load(open(fname, "rb"))

plt.figure()
# Contrast stretching
p1, p2 = np.percentile(img, (0.05, 99))
img_rescale = exposure.rescale_intensity(img, in_range=(p1, p2))
plt.imshow(img_rescale, cmap="gray")

for cnt in contours:
    x = cnt[:, 0][:, 0]
    y = cnt[:, 0][:, 1]
    plt.plot(x, y, "r-", linewidth=0.2)
    # still the last closing line will be missing, get it below
    xlast = [x[-1], x[0]]
    ylast = [y[-1], y[0]]
    plt.plot(xlast, ylast, "r-", linewidth=0.2)

plt.xlim(0, 428)
plt.ylim(0, 684)
plt.tight_layout()
plt.axis("scaled")
plt.axis("off")
path_save = path_img.strip(".tif") + "_ilastik.png"
plt.savefig(path_save, format="png", bbox_inches="tight", dpi=300)
