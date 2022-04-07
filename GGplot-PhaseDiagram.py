import matplotlib.pyplot as plt
import numpy as np

################################
# Imput raw data
Factor_Unit_x = '[FUS], uM'
Factor_Unit_y = '[NaCl], mM'
title = 'No-tag FUS phase diagram 20210826'

LLPSconditions_x = np.array([0.1, 0.5, 1, 2, 5, 0.5, 1, 2, 5, 0.5, 1, 2, 5, 1, 2, 5])
LLPSconditions_y = np.array([50, 50, 50, 50, 50, 100, 100, 100, 100, 250, 250, 250, 250, 500, 500, 500,])

nonLLPSconditions_x = np.array([0.1, 0.1, 0.1, 0.5])
nonLLPSconditions_y = np.array([100, 250, 500, 500])

################################
# plot
plt.figure(figsize=(7, 5), dpi=200)
plt.scatter(LLPSconditions_x, LLPSconditions_y, s=80, facecolors='r', edgecolors='r', label='Two Phases')
plt.scatter(nonLLPSconditions_x, nonLLPSconditions_y, s=80, facecolors='none', edgecolors='b', label='One Phase')
plt.xlabel(Factor_Unit_x, fontsize=12, fontweight="bold")
plt.ylabel(Factor_Unit_y, fontsize=12, fontweight="bold")
plt.legend(loc='lower right', bbox_to_anchor=(1.32, 0),
          ncol=1, fancybox=True, shadow=True)
plt.title(title, fontsize=12, fontweight="bold")
plt.tight_layout()
plt.show()
