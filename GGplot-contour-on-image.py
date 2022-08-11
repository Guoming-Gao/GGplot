from tkinter import filedialog as fd
from os.path import join, dirname, basename
from tifffile import imread
import matplotlib.pyplot as plt
from skimage import exposure
import numpy as np
import cv2

path_mask = "/Volumes/AnalysisGG/PROCESSED_DATA/2022July-RNAinFUS-preliminary/20220712_FLmRNA_10FUS_1Mg_10Dex_noTotR_24C/high_freq_50Hz_both_on_FOV-5-condensate_AveProj_mask.tif"
path_img = "/Volumes/AnalysisGG/PROCESSED_DATA/2022July-RNAinFUS-preliminary/20220712_FLmRNA_10FUS_1Mg_10Dex_noTotR_24C/high_freq_50Hz_both_on_FOV-5-condensate_AveProj.tif"

mask = imread(path_mask)
_, edges = cv2.threshold(mask, 1, 2, 0)
contours, _ = cv2.findContours(edges, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
# print("Total Number of Contours Found: ", str(len(contours)))


plt.figure()
# Contrast stretching
img = imread(path_img)
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
