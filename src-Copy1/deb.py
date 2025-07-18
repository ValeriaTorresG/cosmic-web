import os
import glob
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from itertools import combinations

DATA_DIR    = '../data/results'
CATALOG_DIR = '../data/catalogs_by_type'
ROSETTES    = range(20)
VOID_LIMIT  = -0.9
KNOT_LIMIT  =  0.9

def compute_r(n_data, n_rand):
    r = (n_data - n_rand) / (n_data + n_rand)
    return r

for r in ROSETTES:
    counts_path = os.path.join(DATA_DIR, f'rosette_{r}_counts.txt')
    if not os.path.exists(counts_path):
        print(f"[R{r}] doesnt exist {counts_path}")
        continue
    n_data, n_rand = np.loadtxt(counts_path, unpack=True)
    rvals = compute_r(n_data, n_rand)

    n_void   = np.sum(rvals <= VOID_LIMIT)
    n_sheet  = np.sum((rvals > VOID_LIMIT) & (rvals <= 0))
    n_fila   = np.sum((rvals > 0) & (rvals <= KNOT_LIMIT))
    n_peak   = np.sum(rvals > KNOT_LIMIT)

    print(f"R ",
          f"void={n_void}, sheet={n_sheet}, filament={n_fila}, peak={n_peak}")

    posfile   = os.path.join(DATA_DIR, f'rosette_{r}_pos.txt')
    pairfile  = os.path.join(DATA_DIR, f'rosette_{r}_pairs.txt')
    if not (os.path.exists(posfile) and os.path.exists(pairfile)):
        continue

    pos   = np.loadtxt(posfile)
    pairs = np.loadtxt(pairfile, dtype=int)

    is_void = rvals <= VOID_LIMIT
    void_ids = np.nonzero(is_void)[0]

    mask = np.isin(pairs[:,0], void_ids) & np.isin(pairs[:,1], void_ids)
    subp = pairs[mask]

    all_ids = np.unique(subp.flatten()) if len(subp) else np.array([],dtype=int)
    included = np.zeros_like(all_ids)
    def find_friends(fid):
        group=[]
        idx = np.where(all_ids==fid)[0][0]
        if included[idx]: return []
        included[idx]=1
        group.append(fid)
        friends = list(subp[subp[:,0]==fid,1]) + list(subp[subp[:,1]==fid,0])
        for f2 in friends:
            group += [f2] + find_friends(f2)
        return sorted(set(group))

    groups = []
    for fid in all_ids:
        grp = find_friends(fid)
        if len(grp)>1: #! aqui habia un 4 en el original
            groups.append(grp)
    print(f"       FoF-void-groups (>4 pts): {len(groups)}\n")