import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from itertools import combinations

project_root = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(project_root, 'src')
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

from make_web_catalog import make_catalog

VOID_LIMIT = -0.9
KNOT_LIMIT = 0.9


def save_positions(df, name, out_dir):
    coords = df[['X', 'Y', 'Z']].values
    path = os.path.join(out_dir, f'rosette_{name}_pos.txt')
    np.savetxt(path, coords, fmt='%.6f')
    return path


def save_pairs(df, name, out_dir):
    coords = df[['X', 'Y', 'Z']].values
    tri = Delaunay(coords)
    pairs = set()
    for simplex in tri.simplices:
        for i, j in combinations(simplex, 2):
            a, b = sorted((int(i), int(j)))
            pairs.add((a, b))
    arr = np.array(sorted(pairs), dtype=int)
    path = os.path.join(out_dir, f'rosette_{name}_pairs.txt')
    np.savetxt(path, arr, fmt='%d')
    return path


def compute_counts(df_data, df_rand, seed=42):
    np.random.seed(seed)
    Np = len(df_data)
    sub_rand = df_rand.sample(n=Np, random_state=seed).reset_index(drop=True)
    df_all = pd.concat([
        df_data.assign(RAN=False),
        sub_rand.assign(RAN=True)
    ], ignore_index=True)
    coords = df_all[['X', 'Y', 'Z']].values
    is_data = ~df_all['RAN'].values

    tri = Delaunay(coords)
    neighbors = {i: set() for i in range(len(coords))}
    for simplex in tri.simplices:
        for i, j in combinations(simplex, 2):
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


def save_counts(n_data, n_rand, name, out_dir):
    arr = np.vstack([n_data, n_rand]).T
    path = os.path.join(out_dir, f'rosette_{name}_counts.txt')
    np.savetxt(path, arr, fmt='%d')
    return path


def process_rosette(rosette, base_url, out_dir, void_limit, knot_limit):
    df_data = pd.read_csv(
        os.path.join(base_url, f'ELG_{rosette}_clustering_data.ecsv'),
        comment='#', sep=r'\s+', engine='python')
    df_rand = pd.read_csv(
        os.path.join(base_url, f'ELG_{rosette}_clustering_rand.ecsv'),
        comment='#', sep=r'\s+', engine='python')

    pos_data = save_positions(df_data, f'{rosette}_data', out_dir)
    pos_rand = save_positions(df_rand, f'{rosette}_rand', out_dir)

    pairs_data = save_pairs(df_data, f'{rosette}_data', out_dir)
    pairs_rand = save_pairs(df_rand, f'{rosette}_rand', out_dir)

    n_data, n_rand = compute_counts(df_data, df_rand)
    counts_data = save_counts(n_data, n_rand, f'{rosette}_data', out_dir)

    n_rand_data, n_data_data = compute_counts(df_rand, df_data)
    counts_rand = save_counts(n_rand_data, n_data_data, f'{rosette}_rand', out_dir)

    for wt in ['sheet', 'filament', 'knot']:
        out_cat = os.path.join(out_dir, f'rosette_{rosette}_{wt}.csv')
        make_catalog(
            posfile=pos_data,
            pairfile=pairs_data,
            countfile=counts_data,
            catalogfile=out_cat,
            webtype=wt,
            void_limit=void_limit,
            knot_limit=knot_limit
        )
    out_void = os.path.join(out_dir, f'rosette_{rosette}_voids.csv')
    make_catalog(
        posfile=pos_rand,
        pairfile=pairs_rand,
        countfile=counts_rand,
        catalogfile=out_void,
        webtype='void',
        void_limit=void_limit,
        knot_limit=knot_limit
    )

    print(f'Roette {rosette}: data catalogs + random voids generados.')
    return

def main():
    p = argparse.ArgumentParser(description='Genera web catalogs para DESI (ASTRA)')
    p.add_argument('--base_url', required=True, help='Ruta a data/coord/ (archivos ECSV)')
    p.add_argument('--out_dir', required=True, help='Directorio donde guardar resultados')
    p.add_argument('--void_limit', type=float, default=VOID_LIMIT)
    p.add_argument('--knot_limit', type=float, default=KNOT_LIMIT)
    p.add_argument('--rosette', type=int, required=True, help='ID de rosette (0-19)')
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    process_rosette(
        rosette=args.rosette,
        base_url=args.base_url,
        out_dir=args.out_dir,
        void_limit=args.void_limit,
        knot_limit=args.knot_limit
    )

if __name__ == '__main__':
    main()
