#!/usr/bin/env python
import os
import glob
import argparse
from astropy.table import Table, vstack
import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.coordinates import SkyCoord
import multiprocessing
import time
import gzip
import shutil

# Placeholder para datos combinados de randoms y longitud de data
RAND_TABLE = None
LEN_DATA = None

# Semilla para reproducibilidad global; los procesos se reseedearán en init_worker
np.random.seed(42)

def spherical_to_cartesian_ast(ra, dec, z):
    """
    Convierte coordenadas esféricas (RA, DEC, z) a cartesianas (Mpc) usando astropy SkyCoord.
    """
    ra_vals = ra.data if hasattr(ra, 'data') else ra
    dec_vals = dec.data if hasattr(dec, 'data') else dec
    dc = cosmo.comoving_distance(z)
    coords = SkyCoord(ra=ra_vals*u.deg, dec=dec_vals*u.deg, distance=dc, frame='icrs')
    x, y, zc = coords.cartesian.xyz.to(u.Mpc).value
    return x, y, zc


def init_worker(rand_files, len_data):
    """
    Inicializa cada proceso cargando la tabla combinada de randoms y longitud de datos.
    """
    global RAND_TABLE, LEN_DATA
    # Cargar y apilar todas las tablas random
    tables = [Table.read(rf) for rf in rand_files]
    RAND_TABLE = vstack(tables)
    LEN_DATA = len_data
    # Reseed de numpy por proceso para muestras distintas
    np.random.seed()


def process_random_run(args):
    """
    Para cada corrida: toma LEN_DATA filas aleatorias de RAND_TABLE, convierte y escribe FITS comprimido.
    args: tupla (run_index, out_path)
    """
    run_idx, out_path = args
    # Selección aleatoria sin reemplazo si es posible
    total = len(RAND_TABLE)
    if total >= LEN_DATA:
        idx = np.random.choice(total, size=LEN_DATA, replace=False)
    else:
        idx = np.random.choice(total, size=LEN_DATA, replace=True)
    subset = RAND_TABLE[idx]
    x_r, y_r, z_r = spherical_to_cartesian_ast(subset['RA'], subset['DEC'], subset['Z'])
    zone_rand = Table([
        subset['TARGETID'], subset['Z'], x_r, y_r, z_r
    ], names=('TARGETID', 'Z', 'X', 'Y', 'Zc'))

    # Escribir FITS temporal y comprimir
    tmp_out = out_path.replace('.gz', '')
    zone_rand.write(tmp_out, format='fits', overwrite=True)
    with open(tmp_out, 'rb') as f_in, gzip.open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_out)


def write_zone_data_randoms(base_dir, out_dir, n_runs=100):
    """
    Procesa zonas NGC/SGC: datos + n_runs randoms muestreados de todos los archivos disponibles.
    Imprime tiempo por zona y total.
    """
    os.makedirs(out_dir, exist_ok=True)
    start_all = time.time()

    for hemi in ['NGC', 'SGC']:
        zone_start = time.time()
        # Leer datos
        dat_path = os.path.join(base_dir, f'ELG_LOPnotqso_{hemi}_clustering.dat.fits')
        data = Table.read(dat_path)
        n_data = len(data)

        # Convertir datos a cartesianas y escribir
        x_d, y_d, z_d = spherical_to_cartesian_ast(data['RA'], data['DEC'], data['Z'])
        zone_data = Table([
            data['TARGETID'], data['Z'], x_d, y_d, z_d
        ], names=('TARGETID', 'Z', 'X', 'Y', 'Zc'))
        data_out = os.path.join(out_dir, f'ELG_LOPnotqso_{hemi}_clustering_data.fits.gz')
        tmp_data = data_out.replace('.gz', '')
        zone_data.write(tmp_data, format='fits', overwrite=True)
        with open(tmp_data, 'rb') as f_in, gzip.open(data_out, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(tmp_data)

        # Preparar randoms
        rand_pattern = os.path.join(base_dir, f'ELG_LOPnotqso_{hemi}_*_clustering.ran.fits')
        rand_files = sorted(glob.glob(rand_pattern))
        if not rand_files:
            raise ValueError(f'No se encontraron archivos random para {hemi}')

        # Inicializar pool con todos los randoms y cantidad de datos
        args_list = []
        for run in range(n_runs):
            out_path = os.path.join(out_dir, f'ELG_LOPnotqso_{hemi}_{run}_clustering_rand.fits.gz')
            args_list.append((run, out_path))

        with multiprocessing.Pool(initializer=init_worker, initargs=(rand_files, n_data)) as pool:
            pool.map(process_random_run, args_list)

        elapsed_zone = time.time() - zone_start
        print(f"Zona {hemi} procesada en {elapsed_zone:.2f} segundos")

    total_elapsed = time.time() - start_all
    print(f"Tiempo total de ejecución (ambas zonas): {total_elapsed:.2f} segundos")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Genera datos FITS comprimidos para DR1 (NGC & SGC) muestreando randoms disponibles.'
    )
    parser.add_argument(
        '--base_dir',
        default='/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/raw/',
        help='Directorio de entrada con archivos raw DR1'
    )
    parser.add_argument(
        '--out_dir',
        default='/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/results/',
        help='Directorio de salida para archivos comprimidos'
    )
    parser.add_argument(
        '--n_runs',
        type=int,
        default=100,
        help='Número de iteraciones random por zona'
    )
    args = parser.parse_args()

    write_zone_data_randoms(args.base_dir, args.out_dir, args.n_runs)
