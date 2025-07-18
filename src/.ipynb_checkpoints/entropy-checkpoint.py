import os
import sys
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import time as t
from scipy.spatial import Delaunay
from itertools import combinations

project_root = os.path.dirname(os.path.dirname(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from src.classify_web import astra

NP_RANDOM_SEED = 42
np.random.seed(NP_RANDOM_SEED)

N_ITER = 100
ROSETTE_IDS = range(20)
TYPES = ['void', 'sheet', 'filament', 'knot']


def entropy(p_w):
    norm = -1.0 / np.log2(len(TYPES))
    with np.errstate(divide='ignore', invalid='ignore'):
        s = p_w * np.log2(p_w)
    s = np.nan_to_num(s, nan=0.0)
    return norm * s.sum(axis=1)


def load_data(rosette, base_url):
    data_path = os.path.join(base_url, f'ELG_{rosette}_clustering_data.ecsv')
    rand_path = os.path.join(base_url, f'ELG_{rosette}_clustering_rand.ecsv')
    df_data = pd.read_csv(data_path, comment='#', sep=r'\s+', engine='python')
    df_rand = pd.read_csv(rand_path, comment='#', sep=r'\s+', engine='python')
    return df_data, df_rand


def save_positions(rosette, df_data, out_dir):
    filename = f'rosette_{rosette}_pos.txt'
    path = os.path.join(out_dir, filename)
    coords = df_data[['X', 'Y', 'Z']].values
    np.savetxt(path, coords, fmt='%.6f')
    return path


def compute_counts(df_data, df_rand):
    Np = len(df_data)
    sub_rand = df_rand.sample(n= Np, random_state=NP_RANDOM_SEED).reset_index(drop=True)
    df_all = pd.concat([
        df_data.assign(RAN=False),
        sub_rand.assign(RAN=True)
    ], ignore_index=True)
    coords = df_all[['X','Y','Z']].values
    is_data = ~df_all['RAN'].values

    tri = Delaunay(coords)
    neighbors = {i: set() for i in range(len(coords))}
    for simplex in tri.simplices:
        for i,j in combinations(simplex,2):
            neighbors[i].add(j)
            neighbors[j].add(i)

    n_data = np.zeros(len(coords), dtype=int)
    n_rand = np.zeros(len(coords), dtype=int)
    for i, nbrs in neighbors.items():
        nd = int(np.sum(is_data[list(nbrs)]))
        nr = len(nbrs) - nd
        n_data[i] = nd
        n_rand[i] = nr

    return n_data[:Np], n_rand[:Np]


def save_counts(rosette, n_data, n_rand:, out_dir):
    filename = f'rosette_{rosette}_counts.txt'
    path = os.path.join(out_dir, filename)
    arr = np.vstack([n_data, n_rand]).T
    np.savetxt(path, arr, fmt='%d')
    return path


def save_pairs(rosette, df_data, out_dir):
    coords = df_data[['X','Y','Z']].values
    tri = Delaunay(coords)
    pairs = set()
    for simplex in tri.simplices:
        for i,j in combinations(simplex,2):
            a,b = sorted((int(i),int(j)))
            pairs.add((a,b))
    arr = np.array(sorted(pairs), dtype=int)
    filename = f'rosette_{rosette}_pairs.txt'
    path = os.path.join(out_dir, filename)
    np.savetxt(path, arr, fmt='%d')
    return path


def run_iterations(df_data, df_rand, n_iter, rosette):
    Np = len(df_data)
    counts = np.zeros((Np, len(TYPES)), dtype=int)
    r_list = None

    for it in range(n_iter):
        idx = np.random.choice(len(df_rand), size= Np, replace=False)
        sub_rand = df_rand.iloc[idx].reset_index(drop=True)
        df_typed, _, _, is_real = astra(df_data.copy(), sub_rand.copy())

        labels = df_typed.loc[is_real, 'TYPE'].values
        r_vals = df_typed.loc[is_real, 'r'].values
        if it == 0:
            r_list = r_vals.copy()
        idx_map = [TYPES.index(lbl) for lbl in labels]
        counts[np.arange(Np), idx_map] += 1

    return counts, r_list


def save_results(rosette, df_data, p_w, H, r_list, out_dir):
    final_types = [TYPES[i] for i in np.argmax(p_w, axis=1)]
    df_out = pd.DataFrame({
        'TARGETID': df_data['TARGETID'],
        'r': r_list,
        'type': final_types,
        **{f'p_{t}': p_w[:,i] for i,t in enumerate(TYPES)},
        'H': H
    })
    filename = f'rosette_{rosette}.csv'
    path = os.path.join(out_dir, filename)
    df_out.to_csv(path, index=False)
    return path


def process_rosette(rosette, base_url, out_dir):
    df_data, df_rand = load_data(rosette, base_url)

    pos_path = save_positions(rosette, df_data, out_dir)
    counts_path = save_counts(rosette, *compute_counts(df_data, df_rand), out_dir)
    pairs_path  = save_pairs(rosette, df_data, out_dir)

    counts_iter, r_list = run_iterations(df_data, df_rand, N_ITER, rosette)
    p_w = counts_iter.astype(float)/N_ITER
    H = entropy(p_w)
    results_path = save_results(rosette, df_data, p_w, H, r_list, out_dir)

    print(f'[{rosette}] pos {pos_path}, counts {counts_path}, pairs {pairs_path}')
    return results_path


def main():
    init_t = t.time()
    base_url = os.path.join(project_root, '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/coord/')
    out_dir = os.path.join(project_root, '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/results/')
    os.makedirs(out_dir, exist_ok=True)

    with ProcessPoolExecutor(max_workers = max(1, os.cpu_count()-1)) as exe:
        futures = {exe.submit(process_rosette, r, base_url, out_dir): r for r in ROSETTE_IDS}
        for fut in as_completed(futures):
            r = futures[fut]
            try:
                print(f'Rosette {r} ASTRA {fut.result()}')
            except Exception as e:
                print(f'ERROR rosette {r}: {e}')
    print(f'Elapsed time {init_t - t.time()} s')

if __name__ == '__main__':
    main()