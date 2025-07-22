#!/usr/bin/env python3
import os, sys, time, random, traceback, multiprocessing
import numpy as np, pandas as pd
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
from astropy.table import Table
from scipy.spatial import Delaunay

project_root = os.path.dirname(os.path.dirname(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
from src.classify_web import astra

REGIONS = ['NGC-1', 'NGC-2', 'SGC-3']
N_RUNS = 100
MAX_WORKERS = 50
BASE_DIR = '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/regions/'
OUT_DIR = '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/entropy/'

def seed_everything():
    np.random.seed(42)
    random.seed(42)

def entropy(p_w: np.ndarray) -> np.ndarray:
    norm = -1.0 / np.log2(p_w.shape[1])
    with np.errstate(divide='ignore', invalid='ignore'):
        s = p_w * np.log2(p_w)
    s = np.nan_to_num(s)
    return norm * s.sum(axis=1)

def load_data(region: str, run: int):
    hemi  = region.split('-')[0]
    dpath = os.path.join(BASE_DIR, f'ELG_LOPnotqso_{hemi}_{region}_data.fits.gz')
    rpath = os.path.join(BASE_DIR, f'ELG_LOPnotqso_{hemi}_{region}_{run}_rand.fits.gz')
    tab_d = Table.read(dpath)
    tab_r = Table.read(rpath)
    return pd.DataFrame(tab_d.as_array()), pd.DataFrame(tab_r.as_array())

def process_run(region: str, run: int):
    try:
        df_d, df_r = load_data(region, run)
        coords = np.vstack([df_d[['X','Y','Zc']].values,
                            df_r[['X','Y','Zc']].values])
        is_data = np.hstack([np.ones(len(df_d), bool),
                             np.zeros(len(df_r), bool)])
        pos_file = os.path.join(OUT_DIR, f'{region}_run{run}_pos.txt')
        np.savetxt(pos_file, df_d[['X','Y','Zc']].values, fmt='%.6f')
        tri = Delaunay(coords)
        neighbors = [set() for _ in coords]
        for simplex in tri.simplices:
            for i, j in combinations(simplex, 2):
                neighbors[i].add(j)
                neighbors[j].add(i)
        npts = len(df_d)
        n_data = np.array([sum(is_data[list(nb)]) for nb in neighbors[:npts]], int)
        n_rand = np.array([len(nb) - sum(is_data[list(nb)]) 
                             for nb in neighbors[:npts]], int)
        count_file = os.path.join(OUT_DIR, f'{region}_run{run}_counts.txt')
        np.savetxt(count_file, np.vstack([n_data, n_rand]).T, fmt='%d')

        pairs = set()
        for i, nb in enumerate(neighbors[:npts]):
            for j in nb:
                if j < npts:
                    pairs.add(tuple(sorted((i, j))))
        pair_file = os.path.join(OUT_DIR, f'{region}_run{run}_pairs.txt')
        np.savetxt(pair_file, np.array(sorted(pairs), int), fmt='%d')

        df_typed, _, _, is_real = astra(df_d.copy(), df_r.copy())
        types  = df_typed.loc[is_real, 'TYPE'].values
        r_vals = df_typed.loc[is_real, 'r'].values
        return types, r_vals

    except Exception:
        tb = traceback.format_exc()
        raise RuntimeError(f'[process_run] region={region} run={run} error:\n{tb}')

def process_region(region: str):
    print(f'\nProcessing {region}')
    # prep
    os.makedirs(OUT_DIR, exist_ok=True)
    df_d, _ = load_data(region, 0)
    npts = len(df_d)
    TYPES = ['void','sheet','filament','knot']

    counts    = np.zeros((npts, len(TYPES)), dtype=int)
    r_matrix  = np.zeros((npts, N_RUNS), dtype=float)
    succ_runs = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS,
                             initializer=seed_everything) as exe:
        futures = {exe.submit(process_run, region, run): run
                   for run in range(N_RUNS)}

        for fut in as_completed(futures):
            run = futures[fut]
            try:
                types, r_vals = fut.result()
            except Exception as e:
                print(f'Run {run} failed: {e}', file=sys.stderr)
                continue

            idxs = [TYPES.index(t) for t in types]
            counts[np.arange(npts), idxs] += 1
            r_matrix[:, run] = r_vals
            succ_runs.append(run)

    n_succ = len(succ_runs)
    if n_succ == 0:
        return

    p_w = counts.astype(float) / n_succ
    H = entropy(p_w)
    r_list = r_matrix[:, succ_runs].mean(axis=1)

    df_out = pd.DataFrame({'TARGETID': df_d['TARGETID'],
                           'r': r_list,
                           'type': [TYPES[i] for i in p_w.argmax(axis=1)],
                           'H': H,**{f'p_{t}': p_w[:, j] for j, t in enumerate(TYPES)}})

    out_path = os.path.join(OUT_DIR, f'{region}_entropy.csv')
    df_out.to_csv(out_path, index=False)
    print(f'REgion{region}: {n_succ}/{N_RUNS} runs OK, out in {out_path}')

def main():
    multiprocessing.set_start_method('spawn', force=True)
    seed_everything()
    t0 = time.time()
    for region in REGIONS:
        process_region(region)
    print(f'\nCompleted, elapsed time: {(time.time() - t0):.1f}s!')

if __name__ == '__main__':
    main()