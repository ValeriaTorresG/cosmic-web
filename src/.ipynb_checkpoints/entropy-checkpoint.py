import os
import sys
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import time as t

project_root = os.path.dirname(os.path.dirname(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
from src.classify_web import astra

NP_RANDOM_SEED = 42
N_ITER = 100
ROSETTE_IDS = range(20)
TYPES = ['void','sheet','filament','knot']


def entropy(p_w: np.ndarray):
    norm = -1.0 / np.log2(len(TYPES))
    with np.errstate(divide='ignore', invalid='ignore'):
        s = p_w * np.log2(p_w)
    s = np.nan_to_num(s, nan=0.0)
    return norm * s.sum(axis=1)


def load_data(rosette: int, base_url: str):
    data_path = os.path.join(base_url, f'ELG_{rosette}_clustering_data.ecsv')
    rand_path = os.path.join(base_url, f'ELG_{rosette}_clustering_rand.ecsv')
    df_data = pd.read_csv(data_path, comment='#', sep=r'\s+', engine='python')
    df_rand = pd.read_csv(rand_path, comment='#', sep=r'\s+', engine='python')
    return df_data, df_rand


def run_iterations(df_data, df_rand, n_iter, seed_offset):
    # np.random.seed(NP_RANDOM_SEED + seed_offset)
    Np = len(df_data)
    counts = np.zeros((Np, len(TYPES)), dtype=int)

    for it in range(n_iter):
        # print(f'Iteration {it + 1}/{n_iter} for rosette {seed_offset}')
        idx = np.random.choice(len(df_rand), size= Np, replace=False)
        sub_rand = df_rand.iloc[idx].reset_index(drop=True)

        df_typed, _, _, is_real = astra(df_data.copy(), sub_rand.copy())
        labels = df_typed.loc[is_real, 'TYPE'].values
        r = df_typed.loc[is_real, 'r'].values #!!! TODO falta poner que quede r en el df final
        idx_map = [TYPES.index(lbl) for lbl in labels]

        counts[np.arange(Np), idx_map] += 1

    return counts


def save_results(rosette, df_data, p_w, H, out_dir):
    df_out = pd.DataFrame({'TARGETID': df_data['TARGETID'],
                           **{f'p_{t}': p_w[:, i] for i, t in enumerate(TYPES)},
                           'H': H})
    filename = f'rosette_{rosette}.csv'
    path = os.path.join(out_dir, filename)
    df_out.to_csv(path, index=False)
    return path


def process_rosette(rosette, base_url, out_dir):
    print(f'--------- Processing rosette {rosette} -----------')
    df_data, df_rand = load_data(rosette, base_url)
    counts = run_iterations(df_data, df_rand, N_ITER, seed_offset=rosette)
    p_w = counts.astype(float) / N_ITER
    H = entropy(p_w)
    return save_results(rosette, df_data, p_w, H, out_dir)


def main():
    initial = t.time()
    base_url = os.path.join(project_root, 'data/coord/')
    out_dir  = os.path.join(project_root, 'data/results/')
    os.makedirs(out_dir, exist_ok=True)

    max_workers = max(1, os.cpu_count() - 1)

    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        futures = {exe.submit(process_rosette, r, base_url, out_dir): r
                   for r in ROSETTE_IDS}
        for fut in as_completed(futures):
            r = futures[fut]
            try:
                path = fut.result()
                print(f'Rosette {r} results in {path}')
            except Exception as e:
                print(f'Error {r}: {e}')
    print(f'time {initial-t.time()} sec')


if __name__ == '__main__':
    main()