import numpy as np
import sys
import glob, os
from hyperion_input_seds import get_input_seds
from hyperion.model import ModelOutput
import astropy.constants as constants
import astropy.units as u
import pandas as pd
from tqdm.auto import tqdm
import yt
import fsps



#-------------------------------------------------                                                                                                                                

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

#-------------------------------------------------    

prosp_dir = '/orange/narayanan/s.lower/prospector/mirkwood_comp/simba/z2/snr_10/'
pd_dir = '/orange/narayanan/s.lower/simba/pd_runs/snap160/'
props_dir = '/orange/narayanan/s.lower/simba/ml_files/'


galaxy = sys.argv[1]
sim = sys.argv[2]

print('getting true quantities:')
solar = 0.0196
pdfile = pd_dir+'/grid_physical_properties.160_galaxy'+str(galaxy)+'.npz'
pd_data = np.load(pdfile)
print('    mass and metallicity')
int_mass = np.log10(np.sum(pd_data['grid_star_mass']) / 1.989e33)
int_Z_absolute = np.sum(pd_data['particle_star_metallicity'] * pd_data['particle_star_mass']) / np.sum(pd_data['particle_star_mass'])
int_logzsol = np.log10(int_Z_absolute / solar)

print('    sfr')
print('        loading snap with yt')
ds = yt.load('/orange/narayanan/s.lower/simba/filtered_snapshots/snap160/galaxy_'+str(galaxy)+'.hdf5')
ad = ds.all_data()
star_masses = ad[('PartType4', 'Masses')].in_units('Msun')
scalefactor = ad[("PartType4", "StellarFormationTime")]
star_metal = ad[("PartType4", 'metallicity')]
formation_z = (1.0 / scalefactor) - 1.0
stellar_formation_age = ds.cosmology.t_from_z(formation_z).in_units("Gyr")
simtime = ds.cosmology.t_from_z(ds.current_redshift).in_units("Gyr")
stellar_ages = (simtime - stellar_formation_age).in_units("Gyr")
w50 = np.where(stellar_ages.in_units('Myr').value < 100)[0]
initial_mass = 0.0
print('        initializing sp')
model_sp = fsps.StellarPopulation(zcontinuous=1, sfh=0, add_dust_emission=False, logzsol=0.0, dust2=0.0)

print('        getting initial SP mass')
for star in tqdm(w50):
    current_mass = star_masses[star]
    model_sp.params["logzsol"] = np.log10(star_metal[star] / solar)
    mass_remaining = model_sp.stellar_mass
    initial_mass += current_mass / np.interp(np.log10(stellar_ages[star]*1.e9),model_sp.ssp_ages,mass_remaining)
sfr_100myr = initial_mass/100.e6
print('    dust mass')


if sim == 'simba':
    int_dust_mass = np.sum(pd_data['particle_dustmass'])
else:
    int_dust_mass = 




print('done. saving')
data = {'Galaxy' : galaxy,  'true_mass' : int_mass, 'true_logszol': int_logzsol, 'true_dustmass': int_dust_mass, 'true_sfr': sfr_100myr}


np.savez(prosp_dir+'/galaxy_'+str(galaxy)+'_true_quan.npz', data=data)
