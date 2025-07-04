import os
from astropy.table import Table
import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo


def spherical_to_cartesian(ra, dec, z):
    d_z = cosmo.comoving_distance(z).to(u.Mpc).value  #z to comoving dist in Mpc
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)

    x = d_z * np.cos(dec_rad) * np.cos(ra_rad)
    y = d_z * np.cos(dec_rad) * np.sin(ra_rad)
    z = d_z * np.sin(dec_rad)

    return (x, y, z)


def write_rosette_ecsvs(data, rand_data, rosette_ids, out_dir):
    '''
    Loop over survey data, convert sperical to cartesian coordinates,
    and save data and random catalogs to ECSV.
    '''
    os.makedirs(out_dir, exist_ok=True)

    for rosette in rosette_ids:
        try:
            mask_d = data['ROSETTE_NUMBER'] == rosette
            mask_r = rand_data['ROSETTE_NUMBER'] == rosette
            rosette_dat = data[mask_d]
            rosette_rand = rand_data[mask_r]

            # convert to cart
            x_data, y_data, z_data = spherical_to_cartesian(
                rosette_dat['RA'], rosette_dat['DEC'], rosette_dat['Z']
            )
            x_rand, y_rand, z_rand = spherical_to_cartesian(
                rosette_rand['RA'], rosette_rand['DEC'], rosette_rand['Z']
            )

            # astropy tables to save
            tab_data = Table(
                [rosette_dat['TARGETID'], x_data, y_data, z_data],
                names=('targetid', 'x', 'y', 'z')
            )
            tab_rand = Table(
                [rosette_rand['TARGETID'], x_rand, y_rand, z_rand],
                names=('targetid', 'x', 'y', 'z')
            )

            # write ecsv
            data_out = os.path.join(out_dir, f'ELG_{rosette}_clustering_data.ecsv')
            rand_out = os.path.join(out_dir, f'ELG_{rosette}_clustering_rand.ecsv')

            tab_data.write(data_out)#, overwrite=True)
            tab_rand.write(rand_out)#, overwrite=True)

        except Exception as e:
            # raise ValueError('Error in data')
            print(f'Error with rosette {rosette}: {e}')
            continue