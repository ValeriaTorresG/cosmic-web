import os, time, gzip, shutil, random, sys
from itertools import combinations
import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import Delaunay
from concurrent.futures import ProcessPoolExecutor, as_completed

project_root = os.path.dirname(os.path.dirname(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
from src.classify_web import astra

# Regiones de alta completitud
REGIONS = ['NGC-1', 'NGC-2', 'SGC-3']
# Número de réplicas aleatorias y de iteraciones de clasificación
N_ITER      = 100
MAX_WORKERS = 20

# Directorios (ajusta según tu entorno)
BASE_DIR = '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/regions/'
OUT_DIR  = '/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/entropy/'

# Semillas para reproducibilidad
def seed():
    np.random.seed(42)
    random.seed(42)

# Entropía normalizada
def entropy(p_w: np.ndarray) -> np.ndarray:
    norm = -1.0 / np.log2(len(p_w[0]))
    with np.errstate(divide='ignore', invalid='ignore'):
        s = p_w * np.log2(p_w)
    s = np.nan_to_num(s)
    return norm * s.sum(axis=1)

# Carga las tablas FITS y las convierte en DataFrames
def load_data(region: str, run: int):
    hemi = region.split('-')[0]
    dpath = os.path.join(BASE_DIR, f'ELG_LOPnotqso_{hemi}_{region}_data.fits.gz')
    rpath = os.path.join(BASE_DIR, f'ELG_LOPnotqso_{hemi}_{region}_{run}_rand.fits.gz')
    tab_d = Table.read(dpath);
    tab_r = Table.read(rpath);
    return pd.DataFrame(tab_d.as_array()), pd.DataFrame(tab_r.as_array())

# Procesa una única tarea (región + réplica)
def process_task(args):
    region, run = args
    df_d, df_r = load_data(region, run)
    # Coordenadas y flags
    coords = np.vstack([df_d[['X','Y','Zc']].values,
                        df_r[['X','Y','Zc']].values])
    is_data = np.hstack([np.ones(len(df_d), bool), np.zeros(len(df_r), bool)])

    # 1) Guardar posiciones reales
    pos_file = os.path.join(OUT_DIR, f'{region}_run{run}_pos.txt')
    np.savetxt(pos_file, df_d[['X','Y','Zc']].values, fmt='%.6f')

    # 2) Triangulación única
    tri = Delaunay(coords)
    neighbors = [set() for _ in coords]
    for simplex in tri.simplices:
        for i, j in combinations(simplex, 2):
            neighbors[i].add(j)
            neighbors[j].add(i)

    # 3) Contar vecinos de datos vs. random
    npts = len(df_d)
    n_data = np.array([sum(is_data[list(nb)]) for nb in neighbors[:npts]], int)
    n_rand = np.array([len(nb) - sum(is_data[list(nb)]) for nb in neighbors[:npts]], int)
    count_file = os.path.join(OUT_DIR, f'{region}_run{run}_counts.txt')
    np.savetxt(count_file, np.vstack([n_data, n_rand]).T, fmt='%d')

    # 4) Guardar pares de datos
    pairs = set()
    for i, nb in enumerate(neighbors[:npts]):
        for j in nb:
            if j < npts:
                pairs.add(tuple(sorted((i, j))))
    pair_file = os.path.join(OUT_DIR, f'{region}_run{run}_pairs.txt')
    np.savetxt(pair_file, np.array(sorted(pairs), int), fmt='%d')

    # 5) Iteraciones de clasificación
    counts = np.zeros((npts, len(['void','sheet','filament','knot'])), int)
    for _ in range(N_ITER):
        df_typed, _, _, is_real = astra(df_d.copy(), df_r.copy())
        lbls = df_typed.loc[is_real, 'TYPE'].values
        idxs = [ ['void','sheet','filament','knot'].index(l) for l in lbls ]
        counts[np.arange(npts), idxs] += 1
    p_w = counts.astype(float) / N_ITER
    H   = entropy(p_w)

    # 6) Guardar resultados finales
    r_list = df_typed.loc[is_real, 'r'].values
    df_out = pd.DataFrame({
        'TARGETID': df_d['TARGETID'],
        'r':         r_list,
        'type':      [ ['void','sheet','filament','knot'][i] for i in p_w.argmax(axis=1) ],
        **{f'p_{t}': p_w[:,j] for j,t in enumerate(['void','sheet','filament','knot'])},
        'H': H
    })
    out_file = os.path.join(OUT_DIR, f'{region}_run{run}.csv')
    df_out.to_csv(out_file, index=False)
    return out_file

# Ejecución en paralelo
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    seed()
    tasks = [(r, i) for r in REGIONS for i in range(N_ITER)]
    with ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=seed) as exe:
        futures = {exe.submit(process_task, t): t for t in tasks}
        for fut in as_completed(futures):
            region, run = futures[fut]
            try:
                print(f'[{region} run{run}] ->', fut.result())
            except Exception as e:
                print(f'ERROR [{region} run{run}]:', e)

if __name__ == '__main__':
    main()