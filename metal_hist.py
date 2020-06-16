#!/usr/bin/env python
# coding: utf-8

import yt
import caesar
import numpy as np
import fsps
import sys, os
import tqdm

caesar_snap = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0305_z0.000.hdf5'
#snapshot = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_305.hdf5'

snapshot = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/snap305/galaxy_196.hdf5'

gal_list = np.load('/orange/narayanan/s.lower/simba/snap305_galaxy.npz', allow_pickle=True)

galaxy_num = int(sys.argv[1])

gal = gal_list['galid'][galaxy_num]


obj = caesar.quick_load(caesar_snap)
solar_Z = 0.0196

ds = yt.load(snapshot)
dd = ds.all_data()
#obj.yt_dataset = ds
#dd = obj.yt_dataset.all_data()

scalefactor = dd[("PartType4", "StellarFormationTime")]
star_masses = dd[("PartType4", "Masses")]
star_Z = dd[("PartType4", "metallicity")]
formation_z = (1.0 / scalefactor) - 1.0
yt_cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=0.68)
stellar_formation_age = yt_cosmo.t_from_z(formation_z).in_units("Gyr")
simtime = yt_cosmo.t_from_z(ds.current_redshift).in_units("Gyr")
stellar_ages = (simtime - stellar_formation_age).in_units("Gyr")

#fsps_years = np.linspace(0.01, 14, 200)
#fsps_years = np.arange(0.01, simtime.in_units("Gyr").value, 0.1)
#fsps_years = np.arange(0, simtime.in_units("Gyr").value+0.1, 0.1)
#sfr_bins = fsps_years
#sfr_times = [(j+i)/2 for i, j in zip(sfr_bins[:-1], sfr_bins[1:])] #lookback time
fsps_ssp = fsps.StellarPopulation(sfh=0,
                zcontinuous=1,
                imf_type=2,
                zred=0.0, add_dust_emission=False)

fsps_years = (10**fsps_ssp.ssp_ages) / 1.0e9

#time_list = [14.2 - item for item in sfr_times]
metal_hist = []
#time_list.append(time_list[-1] + 0.1)

for galaxy in [galaxy_num]:
    print('On galaxy: ',galaxy)
    sage = stellar_ages
    smass = star_masses.in_units('Msun')
    stellar_metallicity=star_Z
    for i in range(len(fsps_years)-1):
        age_bin = np.where((sage > (fsps_years[i])) & (sage < (fsps_years[i+1])))[0]
        print(np.shape(age_bin))
        if len(age_bin) == 0:
            #print(sfr_bins[i])
            if fsps_years[i] < 1.:
                metal_hist.append(-6)
            else:
                metal_hist.append(metal_hist[i-1])
        else:
            current_mass = smass[age_bin]
            initial_mass = np.array([])
            Z = np.array([])
            print('on age bin: ',age_bin)
            for star in tqdm.tqdm(age_bin):
                Z = np.append(Z, np.log10(stellar_metallicity[star] / solar_Z))
                fsps_ssp.params["logzsol"] = np.log10(stellar_metallicity[star] / solar_Z)
                mass_remaining = fsps_ssp.stellar_mass
                initial_mass = np.append(initial_mass, smass[star] / np.interp(np.log10(sage[star]*1.e9),fsps_ssp.ssp_ages,mass_remaining))  
            metal_hist.append(np.sum(Z * initial_mass) / np.sum(initial_mass))
            


#data = {'Galaxy': gal, 'Z(t)' : metal_hist, 'Time': time_list}
np.savez('/orange/narayanan/s.lower/simba/metal_hist2/metallicity_snap305_'+str(gal)+'_highres.npz', Galaxy=gal, Z=metal_hist[:-1][::-1], Time=fsps_years[:-2])
#s = pd.DataFrame(data, index=index)
#s.to_pickle('/ufrc/narayanan/s.lower/sfh_snap305_raised_lower_lim.pkl') 


'''c_timelist = list(dict.fromkeys(old_caesar.index.get_level_values('Time [Gyr]')))
c_sfrlist = old_caesar.loc[[galaxy]]['SFR']


ultra_timelist = list(dict.fromkeys(ultra.index.get_level_values('Time [Gyr]')))
ultra_sfrlist = list(dict.fromkeys(ultra.loc[['{0:03}'.format(galaxy)]]['SFR']))
ultra_sfr16 = list(dict.fromkeys(ultra.loc[['{0:03}'.format(galaxy)]]['SFR_16']))
ultra_sfr84 = list(dict.fromkeys(ultra.loc[['{0:03}'.format(galaxy)]]['SFR_84']))


plt.figure(figsize=(10, 8))
plt.rc('axes', linewidth=2)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.tick_params(axis='both', which='minor')
plt.plot(simtime - sfr_bins[:len(sfr)], sfr, lw=1, label='New', zorder=1)
plt.plot(c_timelist, c_sfrlist, lw=1, label='Old', zorder=0)


plt.yscale('log')
plt.ylim([1e-2, 1e3])
plt.xlim([0, 14])
plt.ylabel('SFR', fontsize=23)
plt.xlabel('Time [Gyr]', fontsize=23)
plt.title('SFR: Galaxy '+str(galaxy), fontsize=23)
plt.legend(loc='best', fontsize=15)
plt.savefig('/ufrc/narayanan/s.lower/SEDs/sfhs/sfr_comp'+str(galaxy)+'.png', dpi=300, bbox_inches='tight')'''






