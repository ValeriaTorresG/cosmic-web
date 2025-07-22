import os, glob, argparse, time, gzip, shutil, random
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.coordinates import SkyCoord
from concurrent.futures import ProcessPoolExecutor, as_completed

#takes ~90s
#ngc1 1821322, ngc2 34432, sgc 74059 using 0.9 and not 0.8
REGION_CUTS = {'NGC-1': {'ra_min':110,'ra_max':260,'dec_min':-10,'dec_max':8,'z_min':0.4,'z_max':0.9},
               'NGC-2': {'ra_min':180,'ra_max':260,'dec_min':30,'dec_max':40,'z_min':0.4,'z_max':0.9},
               'SGC-3': {'z_min':0.4,'z_max':0.9}} #orig is 0.4-0.8, but is empty}


def seed_main():
    np.random.seed(42)
    random.seed(42)

def spherical_to_cartesian_ast(ra, dec, z):
    if len(z) == 0: return np.array([]),np.array([]),np.array([])
    ra_vals = ra.data if hasattr(ra,'data') else ra
    dec_vals = dec.data if hasattr(dec,'data') else dec
    dc = cosmo.comoving_distance(z)
    coords = SkyCoord(ra=ra_vals*u.deg, dec=dec_vals*u.deg, distance=dc, frame='icrs')
    x,y,zc = coords.cartesian.xyz.to(u.Mpc).value
    return x,y,zc


def mask_region(table, cuts):
    m = (table['Z'] >= cuts['z_min']) & (table['Z'] < cuts['z_max'])
    if 'ra_min' in cuts:
        m &= (table['RA'] >= cuts['ra_min']) & (table['RA'] < cuts['ra_max'])
        m &= (table['DEC'] >= cuts['dec_min']) & (table['DEC'] < cuts['dec_max'])
    return m


def process_run_region(run_idx, out_path, rand_files, cuts, n_reg):
    file_path = random.choice(rand_files)
    rand_tab = Table.read(file_path)
    m_rand = mask_region(rand_tab,cuts)
    sub_rand = rand_tab[m_rand]
    idx = np.random.choice(len(sub_rand), size=n_reg, replace=len(sub_rand) < n_reg)
    subset = sub_rand[idx]

    x_r, y_r, z_r = spherical_to_cartesian_ast(subset['RA'], subset['DEC'], subset['Z'])
    out_t = Table([subset['TARGETID'], subset['Z'], x_r,y_r,z_r], names=('TARGETID','Z','X','Y','Zc'))

    tmp = out_path.replace('.gz', '')
    out_t.write(tmp, format='fits', overwrite=True)

    with open(tmp,'rb') as fin,gzip.open(out_path,'wb') as fout:shutil.copyfileobj(fin,fout)
    os.remove(tmp)
    return run_idx


def write_data_randoms(base_dir, out_dir, n_runs, max_workers=8):
    os.makedirs(out_dir, exist_ok=True); seed_main()
    t0 = time.time()

    for hemi in ['NGC','SGC']:
        print(f'------------- Zone {hemi}')
        dat = Table.read(os.path.join(base_dir,f'ELG_LOPnotqso_{hemi}_clustering.dat.fits'))
        rand_files = sorted(glob.glob(os.path.join(base_dir,f'ELG_LOPnotqso_{hemi}_*_clustering.ran.fits')))
        print(f'----- {len(dat)} galaxies in total')

        regs = [r for r in REGION_CUTS if r.startswith(hemi)]
        for region in regs:
            cuts = REGION_CUTS[region]
            # print(f' -- subset of {region}: with z in [{cuts["z_min"]},{cuts["z_max"]}]'
            #       + (f', RA [{cuts["ra_min"]},{cuts["ra_max"]}], DEC[{cuts["dec_min"]},{cuts["dec_max"]}]' if 'ra_min' in cuts else ''))
            m = mask_region(dat,cuts)
            dat_reg = dat[m]
            n_reg = len(dat_reg)

            print(f'    {n_reg} galaxies in total')
            if n_reg == 0: print(f'-> No data in region {region}'); continue

            x_d,y_d,z_d = spherical_to_cartesian_ast(dat_reg['RA'], dat_reg['DEC'], dat_reg['Z'])
            data_table = Table([dat_reg['TARGETID'], dat_reg['Z'],x_d,y_d,z_d], names=('TARGETID','Z','X','Y','Zc'))
            fn = os.path.join(out_dir, f'ELG_LOPnotqso_{hemi}_{region}_data.fits.gz')
            tmp = fn.replace('.gz', '')
            data_table.write(tmp,format='fits', overwrite=True)

            with open(tmp,'rb') as fin,gzip.open(fn,'wb') as fout:shutil.copyfileobj(fin,fout) #remove tmp file
            os.remove(tmp)
            print(f'  Generating {n_runs} random samples for region {region}')
            with ProcessPoolExecutor(max_workers=max_workers, initializer=seed_main) as exe:
                futures = {exe.submit(process_run_region,i,os.path.join(out_dir,f'ELG_LOPnotqso_{hemi}_{region}_{i}_rand.fits.gz'),rand_files,cuts,n_reg): i for i in range(n_runs)}
                for fut in as_completed(futures):
                    i = futures[fut]
                    try: fut.result(); print(f'     run {i} finished')
                    except Exception as e: print(f'     ERROR run {i}: {e}')
    print(f'All completed in {time.time()-t0:.1f} s')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--base_dir', default='/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/raw/')
    p.add_argument('--out_dir', default='/pscratch/sd/v/vtorresg/cosmic-web/data/dr1/regions/')
    p.add_argument('--n_runs', type=int, default=100)
    p.add_argument('--max_workers', type=int, default=20)
    args = p.parse_args()
    write_data_randoms(args.base_dir,args.out_dir,args.n_runs,args.max_workers)