import os, glob, argparse
from astropy.table import Table
import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.coordinates import SkyCoord
from concurrent.futures import ProcessPoolExecutor, as_completed
import time, gzip, shutil
import random

#! goes from 55GB to 8.6GB

RAND_FILES = None
LEN_DATA = None

def seed_main():
    np.random.seed(42)
    random.seed(42)


def spherical_to_cartesian_ast(ra, dec, z):
    '''
    Convert spherical coordinates (RA, DEC, Z) to Cartesian coordinates (X, Y, Zc).
    RA and DEC are in degrees, Z is redshift.
    Returns X, Y, Zc in Mpc.
    '''
    ra_vals = ra.data if hasattr(ra, 'data') else ra
    dec_vals = dec.data if hasattr(dec, 'data') else dec
    dc = cosmo.comoving_distance(z)
    coords = SkyCoord(ra = ra_vals*u.deg, dec = dec_vals*u.deg, distance = dc, frame='icrs')
    x, y, zc = coords.cartesian.xyz.to(u.Mpc).value
    return x, y, zc


def init_worker(rand_files, len_data):
    '''
    Initialize worker with random files and length of data.
    '''
    global RAND_FILES, LEN_DATA
    RAND_FILES = rand_files
    LEN_DATA = len_data
    np.random.seed()
    random.seed()


def process_run(run_idx, out_path):
    '''
    Process a single run to generate random samples.
    '''
    file_path = random.choice(RAND_FILES)
    rand_table = Table.read(file_path)
    total = len(rand_table)
    replace = total < LEN_DATA
    idx = np.random.choice(total, size = LEN_DATA, replace = replace)
    subset = rand_table[idx]

    x_r, y_r, z_r = spherical_to_cartesian_ast(subset['RA'], subset['DEC'], subset['Z'])
    out_table = Table([subset['TARGETID'], subset['Z'], x_r, y_r, z_r],
                      names=('TARGETID','Z','X','Y','Zc'))

    tmp = out_path.replace('.gz','')
    out_table.write(tmp, format='fits', overwrite = True)
    with open(tmp,'rb') as fin, gzip.open(out_path,'wb') as fout:
        shutil.copyfileobj(fin, fout)
    os.remove(tmp)
    return run_idx


def write_data_randoms(base_dir, out_dir, n_runs):
    '''
    Write data and random samples to output directory.
    base_dir: Directory containing the input data files.
    out_dir: Directory where the output files will be saved.
    n_runs: Number of random samples to generate.
    '''
    os.makedirs(out_dir, exist_ok=True)
    seed_main()
    start_all = time.time()

    for hemi in ['NGC','SGC']:
        print(f'\nProcessing zone {hemi}')
        zone_start = time.time()

        dat = Table.read(os.path.join(base_dir, f'ELG_LOPnotqso_{hemi}_clustering.dat.fits'))
        n_data = len(dat)
        x_d,y_d,z_d = spherical_to_cartesian_ast(dat['RA'], dat['DEC'], dat['Z'])
        data_table = Table([dat['TARGETID'], dat['Z'], x_d, y_d, z_d],
                           names = ('TARGETID','Z','X','Y','Zc'))

        data_out = os.path.join(out_dir, f'ELG_LOPnotqso_{hemi}_clustering_data.fits.gz')
        tmp_data = data_out.replace('.gz','')
        data_table.write(tmp_data, format = 'fits', overwrite = True)
        with open(tmp_data,'rb') as fin, gzip.open(data_out,'wb') as fout:
            shutil.copyfileobj(fin, fout)
        os.remove(tmp_data)

        rand_pattern = os.path.join(base_dir, f'ELG_LOPnotqso_{hemi}_*_clustering.ran.fits')
        rand_files = sorted(glob.glob(rand_pattern))
        if not rand_files:
            raise RuntimeError(f'No random files for {hemi}')

        print(f'----- Generating {n_runs} ran samples using n = {max(1, os.cpu_count()-1)}')
        init_t = time.time()
        with ProcessPoolExecutor(max_workers = max(1, os.cpu_count()-1),
                                 initializer = init_worker,
                                 initargs = (rand_files, n_data)) as exe:
            futures = {exe.submit(process_run, i,
                                  os.path.join(out_dir, f'ELG_LOPnotqso_{hemi}_{i}_clustering_rand.fits.gz')):
                                      i for i in range(n_runs)}
            for fut in as_completed(futures):
                run_idx = futures[fut]
                try:
                    res = fut.result()
                    print(f'    Run {res} done')
                except Exception as e:
                    print(f'    ERROR run {run_idx}: {e}')

        print(f'Zone {hemi} done in {time.time()-zone_start:.2f} s')

    print(f'\nTotal elapsed: {time.time()-start_all:.2f}s')


if __name__=='__main__':
    parsesr = argparse.ArgumentParser()
    parser.add_argument('--base_dir', default='/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/raw/')
    parser.add_argument('--out_dir', default='/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/results/')
    parser.add_argument('--n_runs', type=int, default=100)
    args = parser.parse_args()
    write_data_randoms(args.base_dir, args.out_dir, args.n_runs)