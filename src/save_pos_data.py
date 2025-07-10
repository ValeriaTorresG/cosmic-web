import os
from astropy.table import Table, vstack
import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.coordinates import SkyCoord

np.random.seed(42)


def spherical_to_cartesian(ra, dec, z):
    r = cosmo.comoving_distance(z).to(u.Mpc).value
    phi = np.radians(ra)
    theta = np.radians(90 - dec)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return (x, y, z)


def spherical_to_cartesian_ast(ra, dec, z):
    ra_deg, dec_deg = ra.data, dec.data
    dc = cosmo.comoving_distance(z)
    coords = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, distance=dc, frame='icrs')
    return coords.cartesian.xyz.to(u.Mpc).value

def write_rosettes(base_dir, out_dir):
    base_dir = '../data/edr/vac/edr/lss/v2.0/LSScats/clustering/'
out_dir = '../data/coord/'
os.makedirs(out_dir, exist_ok=True)

for hemi in ['N', 'S']:

    dat_file = os.path.join(base_dir, f'ELG_{hemi}_clustering.dat.fits')
    data = Table.read(dat_file)
    rosette_ids = np.unique(data['ROSETTE_NUMBER'])

    pattern = os.path.join(base_dir, f'ELG_{hemi}_*_clustering.ran.fits')
    rand_files = sorted(glob.glob(pattern))

    rand_tables = [Table.read(rf) for rf in rand_files]

    for rosette in rosette_ids:

        rosette_dat = data[data['ROSETTE_NUMBER'] == rosette]
        x_d, y_d, z_d = spherical_to_cartesian_ast(
            rosette_dat['RA'], rosette_dat['DEC'], rosette_dat['Z']
        )
        tab_data = Table(
            [rosette_dat['TARGETID'], rosette_dat['Z'], x_d, y_d, z_d],
            names=('TARGETID','z','X','Y','Z')
        )
        tab_data.write(
            os.path.join(out_dir, f'ELG_{hemi}_{rosette}_clustering_data.ecsv'),
            overwrite=True
        )

        partes = [
            rt[rt['ROSETTE_NUMBER'] == rosette]
            for rt in rand_tables
        ]
        partes = [p for p in partes if len(p) > 0]

        if partes:
            all_rand = vstack(partes)
            x_r, y_r, z_r = spherical_to_cartesian_ast(
                all_rand['RA'], all_rand['DEC'], all_rand['Z']
            )
            tab_rand = Table(
                [all_rand['TARGETID'], all_rand['Z'], x_r, y_r, z_r],
                names=('TARGETID','z','X','Y','Z')
            )
        else:
            tab_rand = Table(names=('TARGETID','z','X','Y','Z'))

        tab_rand.write(
            os.path.join(out_dir, f'ELG_{rosette}_clustering_rand.ecsv'),
            overwrite=True
        )

def write_rosette(data, rand_data, rosette_ids, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    for rosette in rosette_ids:
        try:
            mask_d = data['ROSETTE_NUMBER'] == rosette
            mask_r = rand_data['ROSETTE_NUMBER'] == rosette

            rosette_dat = data[mask_d]
            rosette_rand = rand_data[mask_r]

            indices = np.random.choice(len(rosette_rand),
                                       size=len(rosette_dat),
                                       replace=False)
            subset_rand = rosette_rand[indices]

            x_data, y_data, z_data = spherical_to_cartesian(
                rosette_dat['RA'], rosette_dat['DEC'], rosette_dat['Z']
            )
            x_rand, y_rand, z_rand = spherical_to_cartesian(
                subset_rand['RA'], subset_rand['DEC'], subset_rand['Z']
            )

            tab_data = Table(
                [rosette_dat['TARGETID'], x_data, y_data, z_data],
                names=('targetid', 'x', 'y', 'z')
            )
            tab_rand = Table(
                [subset_rand['TARGETID'], x_rand, y_rand, z_rand],
                names=('targetid', 'x', 'y', 'z')
            )

            data_out = os.path.join(out_dir,
                f'ELG_{rosette}_clustering_data.ecsv')
            rand_out = os.path.join(out_dir,
                f'ELG_{rosette}_clustering_rand.ecsv')

            tab_data.write(data_out, overwrite=True)
            tab_rand.write(rand_out, overwrite=True)

        except Exception as e:
            print(f'Error with rosette {rosette}: {e}')
            continue