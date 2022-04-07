import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from progress.bar import IncrementalBar

##############################
# parameters to set

folderpath = '/Users/GGM/Downloads/'
f_rna = 'ACTB  AFTER 2X_CELL-1_right.csv'
f_hops = 'ACTB  AFTER 2X_CELL-1_left.csv'
folderpath_save = folderpath
min_tracklength = 3
r_search_multi = 5    # cut off for RNA search: r_search_multi*estimated R


##############################
# functions
def RNAlocator(df_rna, frame, RNA_IDs):
    # RNA_IDs = np.array([2,7,535], dtype=float)
    df = df_rna[df_rna.FRAME == frame]
    RNA_X = np.array([], dtype=float)
    RNA_Y = np.array([], dtype=float)
    # id = RNA_IDs[0]
    for id in RNA_IDs:
        RNA_X = np.append(RNA_X, df[df.TRACK_ID == id].POSITION_X)
        RNA_Y = np.append(RNA_Y, df[df.TRACK_ID == id].POSITION_Y)
    return RNA_X, RNA_Y


def nearbyRNAs(df_rna, df_hops, r_search_multi, min_tracklength):
    trackids = df_hops.TRACK_ID.unique()
    df = deepcopy(df_hops).reset_index(drop=True)

    lst_frames = np.array([], dtype=float)
    lst_RNAIDs = np.array([], dtype=float)
    lst_RNAX = np.array([], dtype=float)
    lst_RNAY = np.array([], dtype=float)
    lst_condensateIDs = np.array([], dtype=float)
    lst_condensateX = np.array([], dtype=float)
    lst_condensateY = np.array([], dtype=float)
    lst_estR = np.array([], dtype=float)
    lst_distance = np.array([], dtype=float)

    indexes = df.index
    bar = IncrementalBar('Processing...', suffix='%(percent).1f%% - %(eta)ds', max=len(indexes))
    len(indexes)
    # index = indexes[655]
    for index in indexes:
        # fetch a row of a condensate's information
        row = df.iloc[index]
        frame = row.FRAME
        # fetch all RNAs in the same frame as condensate
        RNAinframe_x = df_rna[df_rna.FRAME == frame].POSITION_X.to_numpy()
        RNAinframe_y = df_rna[df_rna.FRAME == frame].POSITION_Y.to_numpy()
        RNAinframe_trackID = df_rna[df_rna.FRAME == frame].TRACK_ID.to_numpy()

        r_search = row.ESTIMATED_DIAMETER * r_search_multi / 2
        distances = np.sqrt((RNAinframe_x-row.POSITION_X)**2 + (RNAinframe_y-row.POSITION_Y)**2)

        N = distances[distances < r_search].size

        if N <= min_tracklength:
            continue

        # Saving distance
        lst_distance = np.append(lst_distance, distances[distances < r_search])
        # Saving RNA trackID
        RNA_IDs = RNAinframe_trackID[distances < r_search]
        lst_RNAIDs = np.append(lst_RNAIDs, RNA_IDs)
        # Saving RNA x and y
        RNA_X, RNA_Y = RNAlocator(df_rna, frame, RNA_IDs)
        lst_RNAX = np.append(lst_RNAX, RNA_X)
        lst_RNAY = np.append(lst_RNAY, RNA_Y)
        # Saving time / frame
        lst_frames = np.append(lst_frames, np.repeat(frame, N))
        # Saving condensate trackID
        lst_condensateIDs = np.append(lst_condensateIDs, np.repeat(row.TRACK_ID, N))
        # Saving condensate x and y
        lst_condensateX = np.append(lst_condensateX, np.repeat(row.POSITION_X, N))
        lst_condensateY = np.append(lst_condensateY, np.repeat(row.POSITION_Y, N))
        # Saving condensate estimated R
        lst_estR = np.append(lst_estR, np.repeat(row.ESTIMATED_DIAMETER/2, N))
        bar.next()

    bar.finish()

    df_out = pd.DataFrame(dtype=float)
    df_out['t'] = lst_frames
    df_out['RNA_trackID'] = lst_RNAIDs
    df_out['RNA_X'] = lst_RNAX
    df_out['RNA_Y'] = lst_RNAY
    df_out['condensate_trackID'] = lst_condensateIDs
    df_out['condensate_X'] = lst_condensateX
    df_out['condensate_Y'] = lst_condensateY
    df_out['est_R'] = lst_estR
    df_out['distance'] = lst_distance

    df_out = df_out.sort_values(by=['RNA_trackID', 't', 'condensate_trackID'])

    return df_out


# main:
os.chdir(folderpath)
df_rna = pd.read_csv(f_rna)
df_hops = pd.read_csv(f_hops)

# Get RNA and HOPS tracks under the current fname title

df_rna = df_rna[['TRACK_ID', 'FRAME', 'POSITION_X', 'POSITION_Y']]
df_rna = df_rna[df_rna.TRACK_ID != 'None']
df_rna = df_rna.astype(float).sort_values(by=['TRACK_ID', 'FRAME'])
df_rna.to_csv('test-rna.csv')

df_hops = df_hops[['TRACK_ID', 'FRAME', 'POSITION_X', 'POSITION_Y', 'ESTIMATED_DIAMETER']]
df_hops = df_hops[df_hops.TRACK_ID != 'None']
df_hops = df_hops.astype(float).sort_values(by=['TRACK_ID', 'FRAME'])
df_hops.to_csv('test-hops.csv')

# df_out = nearbyRNAs(df_rna, df_hops, r_search_multi, min_tracklength)
# os.chdir(folderpath_save)
# fname_save = fname_title + '_nearbyRNAs.csv'
# df_out.to_csv(fname_save, index=False)
