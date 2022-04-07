import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
from scipy.optimize import curve_fit
import seaborn as sns
sns.set(color_codes=True, style='white')

###############################################
# Loading Data
fpath1 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210612-dwelltimestoplot.csv'
fpath2 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20210909-dwelltimestoplot.csv'
fpath3 = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD/20211007-dwelltimestoplot.csv'
folderpath_save = '/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/From_AnalysisSSD'

df1 = pd.read_csv(fpath1)
df2 = pd.read_csv(fpath2)
df3 = pd.read_csv(fpath3)

# Omit Xin's L941 data:
df3 = df3[df3.RNAtype != 'L941']

df_all = pd.concat([df1,df2,df3])
df_all = df_all[df_all.HOPScondition == '300 mM Na+']
df_all.dwelltimes = df_all.dwelltimes.to_numpy(dtype=float)*0.1
df_all = df_all[(df_all.dwelltimes>0) & (df_all.dwelltimes<=20)]

def CDF(t,k):
    return 1-np.exp(-k*t)

def calc_R2(ydata, yfit):
    # residual sum of squares (ss_tot)
    residuals = ydata- yfit
    ss_res = np.sum(residuals**2)
    # total sum of squares (ss_tot)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    # r_squared-value
    r_squared = 1 - (ss_res / ss_tot)

    return r_squared


lst_RNA = ['miR-21ds','miR-21gs','THOR','THOR-d','L941','ActB',"SOX2"]
lst_halflife = np.array([],dtype=float)
lst_R2 = np.array([],dtype=float)

os.chdir(folderpath_save)
i = 3
for i in range(len(lst_RNA)):
    plt.figure(figsize=(9, 3), dpi=300)
    data = df_all[df_all.RNAtype==lst_RNA[i]].dwelltimes

    # for better CDF fitting, filter out the most static portion
    data = data[data<19.9]

    hist, bin_edges, _ = plt.hist(data, bins=100, range=(0,20), density=True, histtype='stepfilled', cumulative=1, color=sns.color_palette()[i], label=lst_RNA[i], alpha=0.5)

    bin_centers = bin_edges[:-1]+0.1
    k_fit, _ = curve_fit(CDF, bin_centers, hist)
    xdata = bin_centers
    ydata = hist
    yfit = CDF(xdata, k_fit[0])
    R2 = calc_R2(ydata, yfit)

    plt.plot(np.arange(0,20,0.05), 1-np.exp(-k_fit[0]*np.arange(0,20,0.05)), color=sns.color_palette()[i], linestyle='-')
    # Assuming first order reaction, also note in python np.log is ln
    halflife = np.log(2) / k_fit[0]

    plt.legend(loc='lower right')
    plt.xlabel('Dwell Time (s)')
    plt.ylabel('CDF')
    plt.xlim(0,20)
    plt.tight_layout()

    comment = 'Dwelling Half-Life: ' + str(round(halflife,2)) + ' s'
    plt.text(19.8, 0.8, comment, horizontalalignment='right', size=20, fontweight='normal')

    comment = 'R^2 = ' + str(round(R2,3))
    plt.text(19.8, 0.6, comment, horizontalalignment='right', size=25, fontweight='bold')

    fname_save = 'Dwell Time Distribution-' + lst_RNA[i] + '-CDF-single_fit.png'
    plt.savefig(fname_save, format='png')
    plt.close()

    lst_halflife = np.append(lst_halflife, halflife)
    lst_R2 = np.append(lst_R2, R2)

df_save = pd.DataFrame(dtype=object)
df_save['RNA'] = lst_RNA
df_save['halflife'] = lst_halflife
df_save['R2'] = lst_R2

df_save.to_csv('Dwell Time fitting results_single.csv', index=False)
