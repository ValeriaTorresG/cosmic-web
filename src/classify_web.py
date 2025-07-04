import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from itertools import combinations

def compute_r(df):
    '''
    Compute r for each point in the dataset using
    r = (n_data - n_rand) / (n_data + n_rand)
    '''
    coords  = df[['X', 'Y', 'Z']].values
    is_data = ~df['RAN'].values

    tri = Delaunay(coords)

    #! adjacency list for neighbors
    neighbors = {i: set() for i in range(len(coords))}
    for simplex in tri.simplices:
        for i, j in combinations(simplex, 2):
            neighbors[i].add(j)
            neighbors[j].add(i)

    r = np.zeros(len(coords), dtype=float)
    for i, nbrs in neighbors.items():
        n_data = int(np.sum(is_data[list(nbrs)]))
        n_rand = len(nbrs) - n_data
        if (n_data + n_rand) > 0:
            r[i] = (n_data - n_rand) / (n_data + n_rand)
        else:
            r[i] = 0.0

    out = df.copy()
    out['r'] = r
    return out

def classify_r(df):
    '''
    Classify points based on the computed r value
    '''
    r = df['r'].values
    conds = [(r >= -1.0) & (r <= -0.9),
             (r > -0.9) & (r <= 0.0),
             (r > 0.0) & (r <= 0.9),
             (r > 0.9) & (r <= 1.0),]
    choices = ['void', 'sheet', 'filament', 'knot']
    df = df.copy()
    df['TYPE'] = np.select(conds, choices, default='error')
    return df