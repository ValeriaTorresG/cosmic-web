import os, sys, gzip, shutil
import numpy as np
import pandas as pd
import time as t
from astropy.table import Table
from scipy.spatial import Delaunay
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

mp.set_start_method('spawn', force=True)

project_root = os.path.dirname(os.path.dirname(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
from src.classify_web import astra

NP_RANDOM_SEED = 42
N_ITER = 100
HEMIS = ['NGC', 'SGC']
TYPES = ['void', 'sheet', 'filament', 'knot']

np.random.seed(NP_RANDOM_SEED)

def entropy(p_w):
    norm = -1.0 / np.log2(len(TYPES))
    with np.errstate(divide='ignore', invalid='ignore'):
        s = p_w * np.log2(p_w)
    s = np.nan_to_num(s, nan=0.0)
    return norm * s.sum(axis=1)


def load_data(hemi, run, base_url):
    data_fits = os.path.join(base_url, f'ELG_LOPnotqso_{hemi}_clustering_data.fits.gz')
    rand_fits = os.path.join(base_url, f'ELG_LOPnotqso_{hemi}_{run}_clustering_rand.fits.gz')
    tbl_data = Table.read(data_fits)
    tbl_rand = Table.read(rand_fits)
    return tbl_data.to_pandas(), tbl_rand.to_pandas()


def save_positions(hemi, run, df, out_dir):
    path = os.path.join(out_dir, f'{hemi}_run{run}_pos.txt')
    np.savetxt(path, df[['X','Y','Z']].values, fmt='%.6f')
    return path


def compute_counts(df, dr):
    Np = len(df)
    sr = dr.sample(n=NP_RANDOM_SEED if False else Np, replace=False)
    all_df = pd.concat([df.assign(RAN=False), sr.assign(RAN=True)], ignore_index=True)
    coords = all_df[['X','Y','Z']].values
    is_data = ~all_df['RAN'].values
    tri = Delaunay(coords)
    neighbors = {i: set() for i in range(len(coords))}
    for simplex in tri.simplices:
        for i, j in combinations(simplex, 2):
            neighbors[i].add(j)
            neighbors[j].add(i)
    n_data = np.array([sum(is_data[list(neighbors[i])]) for i in range(len(coords))])
    n_rand = np.array([len(neighbors[i]) - n_data[i] for i in range(len(coords))])
    return n_data[:Np], n_rand[:Np]


def save_counts(hemi, run, n_data, n_rand, out_dir):
    path = os.path.join(out_dir, f'{hemi}_run{run}_counts.txt')
    np.savetxt(path, np.vstack([n_data, n_rand]).T, fmt='%d')
    return path


def save_pairs(hemi, run, df, out_dir):
    coords = df[['X','Y','Z']].values
    tri = Delaunay(coords)
    pairs = set()
    for simplex in tri.simplices:
        for i, j in combinations(simplex, 2):
            pairs.add(tuple(sorted((i, j))))
    path = os.path.join(out_dir, f'{hemi}_run{run}_pairs.txt')
    np.savetxt(path, np.array(sorted(pairs)), fmt='%d')
    return path


def run_iterations(df: pd.DataFrame, dr):
    Np = len(df)
    counts = np.zeros((Np, len(TYPES)), int)
    r_list = None
    for it in range(N_ITER):
        print(N_ITER)
        idx = np.random.choice(len(dr), Np, replace=False)
        sub = dr.iloc[idx].reset_index(drop=True)
        df_typed, _, _, is_real = astra(df.copy(), sub.copy())
        labels = df_typed.loc[is_real, 'TYPE'].values
        r_vals = df_typed.loc[is_real, 'r'].values
        if r_list is None:
            r_list = r_vals.copy()
        for i, lbl in enumerate(labels):
            counts[i, TYPES.index(lbl)] += 1
    return counts, r_list


def save_results(hemi: str, run: int, df: pd.DataFrame, p_w, H, r_list, out_dir: str):
    final_types = [TYPES[i] for i in np.argmax(p_w, axis=1)]
    df_out = pd.DataFrame({'TARGETID': df['TARGETID'],
                           'r': r_list, 'type': final_types,
                           **{f'p_{t}': p_w[:, i] for i, t in enumerate(TYPES)},
                           'H': H})
    path = os.path.join(out_dir, f'{hemi}_run{run}.csv')
    df_out.to_csv(path, index=False)
    return path


def process_zone_run(hemi: str, run: int, base_url: str, out_dir: str):
    df_data, df_rand = load_data(hemi, run, base_url)
    pos = save_positions(hemi, run, df_data, out_dir)
    cnt = save_counts(hemi, run, *compute_counts(df_data, df_rand), out_dir)
    pr = save_pairs(hemi, run, df_data, out_dir)
    counts_iter, r_list = run_iterations(df_data, df_rand)
    p_w = counts_iter / N_ITER
    H = entropy(p_w)
    res = save_results(hemi, run, df_data, p_w, H, r_list, out_dir)
    return f"pos:{pos},counts:{cnt},pairs:{pr},res:{res}"


def main():
    start = t.time()
    base_dir = '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/results/'
    out_dir = '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/entropy/'
    os.makedirs(out_dir, exist_ok=True)

    for hemi in HEMIS:
        tasks = [(hemi, run, base_dir, out_dir) for run in range(N_ITER)]
        max_workers = 50#min(4, os.cpu_count() - 1)
        print(f"Processing zone {hemi} with {len(tasks)} tasks on {max_workers} workers...")
        with ProcessPoolExecutor(max_workers=max_workers, mp_context=mp.get_context('spawn')) as executor:
            futures = {executor.submit(process_zone_run, *args): args[1] for args in tasks}
            for fut in as_completed(futures):
                run = futures[fut]
                try:
                    outp = fut.result()
                    print(f'[{hemi} run{run}] done in {outp}')
                except Exception as e:
                    print(f'ERROR [{hemi} run{run}]: {e}')

    print(f"Total time: {t.time() - start:.2f}s")

if __name__ == '__main__':
    main()