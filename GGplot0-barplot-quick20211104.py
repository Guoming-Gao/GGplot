import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
sns.set(color_codes=True, style='white')


lst_RNA = ['miR-21ds', 'miR-21gs', 'THOR', 'THOR-d', 'L941', 'ActB', 'SOX2']
static_percent = [0,0,2.16,6.33,6.34,0.06,3.02]
folderpath_save = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD'

os.chdir(folderpath_save)
plt.figure(figsize=(4, 5), dpi=300)
sns.barplot(x=lst_RNA, y=static_percent)
plt.ylabel('Static Dwell State Percentage (%)')
plt.xticks(rotation=45)
plt.xlabel('')
plt.tight_layout()
plt.savefig('Static Dwell State Percentage.png',format='png')
