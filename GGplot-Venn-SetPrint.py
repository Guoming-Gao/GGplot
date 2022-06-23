from os.path import join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles


path_dict = "/Users/GGM/Documents/GGscripts/HOPS-simulation-project/"
df_plot = pd.read_excel(join(path_dict, "result-7-with EnigmRBP-202101updated.xlsx"))
df_plot["Ameya_top_Nmer"] = df_plot["Ameya_top_Nmer"].astype(str)  # debug, IMPORTANT
# print(df_plot.columns)

# Extract subsets as dataframes
df_HOPS = df_plot[df_plot["HOPS?(y/n)"] == "y"]
df_nonHOPS = df_plot[df_plot["HOPS?(y/n)"] == "n"]


df_multimer_Nincluded = df_plot[
    (df_plot["Ameya_top_Nmer"] != "1") & (df_plot["Ameya_top_Nmer"] != "2")
]
df_multimer_withoutN = df_plot[
    (df_plot["Ameya_top_Nmer"] != "1")
    & (df_plot["Ameya_top_Nmer"] != "2")
    & (df_plot["Ameya_top_Nmer"] != "n")
]


df_enigmRBP = df_plot[df_plot["EnigmRBP"] == "y"]
df_RBDrecognized = df_plot[df_plot["RBD_recognized"] == "y"]
df_PDBavailable = df_plot[df_plot["Representative PDB"] != "None"]
df_wholeinteractome = df_plot[df_plot["All_mRNA_interactome"] == "y"]
df_3mer = df_plot[df_plot["Ameya_top_Nmer"] == "3"]


# Convert dataframes into sets
set_HOPS = set(df_HOPS["Entry"])
set_nonHOPS = set(df_nonHOPS["Entry"])
set_multimer_Nincluded = set(df_multimer_Nincluded["Entry"])
set_multimer_withoutN = set(df_multimer_withoutN["Entry"])
# set_enigmRBP = set(df_enigmRBP['Entry'])
set_RBDrecognized = set(df_RBDrecognized["Entry"])
set_PDBavailable = set(df_PDBavailable["Entry"])
set_wholeinteractome = set(df_wholeinteractome["Entry"])
set_enigmRBP = set_wholeinteractome.difference(set_RBDrecognized)
set_3mer = set(df_3mer["Entry"])
df_plot.keys()

# Print out lists of proteins


def print_proteins(df_plot, set_toprint):
    print("Total:", len(set_toprint))
    for key in list(set_toprint):
        """
        if df_plot.loc[key, 'Representative PDB'] != 'None':
            print(df_plot.loc[key, 'SearchKey'],
                  '\tUniProt ID:', key,
                  '\tAmeya_top_Nmer:', df_plot.loc[key, 'Ameya_top_Nmer'])
        else:
            print(df_plot.loc[key, 'SearchKey'],
                  '\tUniProt ID:', key,
                  '\tPDB Not Available')
        """
        print(
            df_plot.loc[key, "SearchKey"],
            "\tUniProt ID:",
            key,
            "\tAmeya_top_Nmer:",
            df_plot.loc[key, "Ameya_top_Nmer"],
        )


df_plot = df_plot.set_index("Entry")
# set_toprint = set_HOPS.intersection(set_multimer)
set_toprint = set_nonHOPS.intersection(set_multimer_Nincluded)
set_toprint = set_toprint.intersection(set_enigmRBP)
# set_toprint = set_HOPS.intersection(set_multimer_withoutN)
# print('\nThe set to print:')
print_proteins(df_plot, set_toprint)

plt.figure(dpi=300, linewidth=10)
venn2([set_HOPS, set_multimer_Nincluded], set_labels=("HOPS", "Multimer"))
plt.title("\u03A9: Screening of 100 human proteins", weight="bold", name="Helvetica")
plt.tight_layout()
plt.savefig(
    "/Users/GGM/Documents/Graduate_Work/PhD-Year4/RNA Society 2022/Venn.svg",
    format="svg",
)


"""
def print_proteins(df_plot, set_toprint):
    print('Total:', len(set_toprint))
    for key in list(set_toprint):
        if df_plot.loc[key, 'Representative PDB'] != 'None':
            print('Gene Name:', df_plot.loc[key, 'SearchKey'],
                  '\tUniProt ID:', key,
                  '\tRepresentative PDB ID:', df_plot.loc[key, 'Representative PDB'])
        else:
            print('Gene Name:', df_plot.loc[key, 'SearchKey'],
                  '\tUniProt ID:', key,
                  '\tRepresentative PDB Not Available')

# non-multimer HOPS proteins with neither enigmRBP nor known RBD
set_nonMulHOPS = set_HOPS.difference(set_multimer)
set_toprint = set_nonMulHOPS.difference(set_enigmRBP)
set_toprint = set_toprint.difference(set_RBDrecognized)
print('\nMonomer or Dimer HOPS proteins with neither enigmRBP nor known RBD:')
print_proteins(df_plot, set_toprint)

# non-multimer HOPS proteins with enigmRBP
set_toprint = set_nonMulHOPS.intersection(set_enigmRBP)
print('\n \nMonomer or Dimer HOPS proteins that are enigmRBP:')
print_proteins(df_plot, set_toprint)

# non-multimer HOPS proteins with known RBD
set_toprint = set_nonMulHOPS.intersection(set_RBDrecognized)
print('\n \nMonomer or Dimer HOPS proteins with known RBD:')
print_proteins(df_plot, set_toprint)

# non-multimer non-HOPS proteins with neither enigmRBP nor known RBD
set_nonMul_nonHOPS = set_nonHOPS.difference(set_multimer)
set_toprint = set_nonMul_nonHOPS.difference(set_enigmRBP)
set_toprint = set_toprint.difference(set_RBDrecognized)
print('\nMonomer or Dimer non-HOPS proteins with neither enigmRBP nor known RBD:')
print_proteins(df_plot, set_toprint)

# non-multimer non-HOPS proteins with enigmRBP
set_toprint = set_nonMul_nonHOPS.intersection(set_enigmRBP)
print('\n \nMonomer or Dimer non-HOPS proteins that are enigmRBP:')
print_proteins(df_plot, set_toprint)

# non-multimer non-HOPS proteins with known RBD
set_toprint = set_nonMul_nonHOPS.intersection(set_RBDrecognized)
print('\n \nMonomer or Dimer non-HOPS proteins with known RBD:')
print_proteins(df_plot, set_toprint)

"""


"""
# Convert list to dict, and plot Venn with venn (UNUSED)
from venn import venn
dict_venn = {
    'HOPS': set_HOPS,
    'Multimer': set_multimer,
    'EnigmRBP': set_enigmRBP,
    'RBD recognized': set_RBDrecognized
}
plt.figure(num=1, figsize=(3, 4), dpi=200)
venn(dict_venn, cmap='Set1')
"""

"""
# plot Venn diagram with venn3 from matplotlib
plt.figure(num=2, figsize=(3, 3), dpi=200)
venn3([set_HOPS, set_multimer, set_enigmRBP],
      set_labels=('HOPS', 'Multimer', 'EnigmRBP'))
plt.figure(num=3, figsize=(3, 3), dpi=200)
venn3([set_HOPS, set_multimer, set_RBDrecognized],
      set_labels=('HOPS', 'Multimer', 'RBD recognized'))

# plot multimer vs mRNA interactome
plt.figure(num=1, figsize=(3, 3), dpi=200)
venn2([set_HOPS, set_multimer],
      set_labels=('HOPS', 'Multimer'))
plt.figure(num=2, figsize=(3, 3), dpi=200)
venn2([set_HOPS, set_wholeinteractome],
      set_labels=('HOPS', 'mRNA\nInteractome'))
plt.figure(num=3, figsize=(3, 3), dpi=200)
venn3([set_HOPS, set_multimer, set_wholeinteractome],
      set_labels=('HOPS', 'Multimer', 'mRNA Interactome'))

plt.show()
"""
