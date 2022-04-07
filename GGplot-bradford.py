import scipy.stats as stats
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd

##############################
# parameters
folderpath = "/Users/GGM/Documents/Graduate_Work/Nils_Walter_Lab/[wLiuhan] FUS FL purification-in vitro SM/protein quantification/"
fname = "20210821-final-formatted.csv"
# FUS 1 mg/ml = 18.72 uM (53.425 kDa)
sample_uM_per_mgml = 18.72
RHKrate_BSA = 0.17
RHKrate_sample = 0.103  # FUS, no tag
convert_ratio = RHKrate_BSA/RHKrate_sample

# Explanation:
# [sample_read] is the sample concentration **assuming** it is BSA.
# [sample_true] * RHKrate_sample = [sample_read] * RHKrate_BSA
# Therefore,
# [sample_true]
# = [sample_read] * (RHKrate_BSA/RHKrate_sample)
# = [sample_read] * convert_ratio


##############################
# Main
fpath = folderpath + fname

df = pd.read_csv(fpath)
df = df.set_index('title')

df['A595-ave'] = df.iloc[:, 0:3].mean(axis=1)
df['A595-std'] = df.iloc[:, 0:3].std(axis=1)
stdcurve_ave = np.array(df['A595-ave'][0:6])
stdcurve_std = np.array(df['A595-std'][0:6])
stdcurve_conc = np.array(df['conc, mg/ml'][0:6])
slope, intercept, r_value, p_value, std_err = stats.linregress(stdcurve_conc, stdcurve_ave)
R2 = r_value**2
line1 = 'slope = ' + str(slope.round(3))
line2 = 'intercept = ' + str(intercept.round(3))
line3 = 'R^2 = ' + str(R2.round(3))
textbox = line1 + '\n' + line2 + '\n' + line3

plt.errorbar(x=stdcurve_conc, y=stdcurve_ave, yerr=stdcurve_std, fmt='b.', capsize=5)
plt.plot(stdcurve_conc, slope*stdcurve_conc+intercept, 'r')
plt.xlabel('Concentration, mg/mL')
plt.ylabel('A595, a.u.')
plt.title(fname[0:-4])
plt.tight_layout()
plt.annotate(textbox, xy=(0.05, 0.8), xycoords='axes fraction', fontsize=13)

os.chdir(folderpath)
fname_fig = fname[0:-4] + '.png'
plt.savefig(fname_fig, format='png')


##############################
# calculate samples and save results
def sampleprint(file, sample, A595, sample_read, sample_true, sample_true_mol):
    file.write('Sample Name: ' + str(sample) + '\n')
    file.write('Average A595 (n=3): ' + str(A595.round(3)) + '\n')
    file.write('Read concentration (mg/mL): ' + str(sample_read.round(3)) + '\n')
    file.write('True concentration (mg/mL): ' + str(sample_true.round(3)) + '\n')
    file.write('Molar concentration (uM): ' + str(sample_true_mol.round(3)) + '\n')
    file.write('\n\n')


fname_txt = fname[0:-4] + '_results' + '.txt'
file = open(fname_txt, "w")

samples = df.index[-4:]
for sample in samples:
    A595 = df.loc[sample, 'A595-ave']
    # A = k * Conc + b
    # Conc = (A - b) / k = (A595-ave - intercept) / slope
    sample_read = (A595 - intercept) / slope
    sample_true = sample_read * convert_ratio
    sample_true_mol = sample_true * sample_uM_per_mgml
    sampleprint(file, sample, A595, sample_read, sample_true, sample_true_mol)

file.close()
