import numpy as np
from scipy.spatial import Delaunay
import pandas as pd
from itertools import combinations


def compute_r(df):
    coords = df[['X', 'Y', 'Z']].values
    is_data = ~df['RAN'].values

    tri = Delaunay(coords)
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
            raise ValueError(f'No neighbors for point {i} in the triangulation.')

    out = df.copy()
    out['r'] = r
    return out


def classify_r(df):
    r = df['r'].values
    conds = [(r >= -1.0) & (r <= -0.9),
             (r > -0.9) & (r <= 0.0),
             (r > 0.0) & (r <= 0.9),
             (r > 0.9) & (r <= 1.0),]
    choices = ['void', 'sheet', 'filament', 'knot']
    df = df.copy()
    df['TYPE'] = np.select(conds, choices, default='error')
    return df


def astra(data, rand):
    data['RAN'] = False
    rand['RAN'] = True
    df = pd.concat([data, rand], ignore_index=True)

    df_r = compute_r(df)
    df_typed = classify_r(df_r)
    coords2d = df_typed[['X','Y','Z']].values
    tri = Delaunay(coords2d)

    is_real = ~df_typed['RAN'].values
    return df_typed, coords2d, tri, is_real