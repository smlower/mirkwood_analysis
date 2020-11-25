from corner import quantile
import prospect.io.read_results as pread
from prospect.models.transforms import zfrac_to_sfrac, logsfr_ratios_to_sfrs, zfrac_to_sfr, tburst_from_fage, zfrac_to_masses, logsfr_ratios_to_masses
import numpy as np
import sys
import glob, os
import pandas as pd

def nonpara_massweighted_age(sfr_in_each_bin, time_bins):
    top = 0.0
    bottom = 0.0
    time = (10**time_bins) / 1e9

    for bin_ in range(len(sfr_in_each_bin)):
        top += np.abs(time[bin_] * sfr_in_each_bin[bin_])
        bottom += np.abs(sfr_in_each_bin)
    return top / bottom

prosp_dir = '/orange/narayanan/s.lower/prospector/mirkwood_comp/eagle/'
#pd_dir = '/orange/narayanan/s.lower/simba/pd_runs/snap305/'
pd_dir = '/orange/narayanan/s.lower/eagle/pd_runs/snap28/'
galaxy_idx = sys.argv[1]

b = pd.read_pickle('/orange/narayanan/s.lower/eagle/ml_files/eagle_ml_SEDs_z0.0.pkl')

galaxy = int(b['ID'].loc[[int(galaxy_idx)]])

#galaxy = np.load('/orange/narayanan/s.lower/prospector/simba_runs/simba_galaxy_SFRcut.npz')['ID'][int(galaxy_idx)]



#appending path of prosp run script to import model and sps
sys.path.append(prosp_dir)
from run_prosp import build_model, build_sps



sfr_50 = [] 
sfr_16 = []
sfr_84 = []

print('now reading files')
galaxy_num = "{:03d}".format(galaxy)
infile = prosp_dir+'/galaxy_'+str(int(galaxy_idx))+'_SNR_10.h5'

print('running eagle galaxy '+str(galaxy)+', prospector galaxy '+str(galaxy_idx))


for prosp_output in glob.glob(infile):
    print(prosp_output)
    res, obs, mod = pread.results_from(prosp_output)
pdfile = pd_dir+'/grid_physical_properties.028_galaxy'+str(galaxy)+'.npz'
pd_data = np.load(pdfile)
int_mass = np.log10(np.sum(pd_data['grid_star_mass']) / 1.989e33)

print('sps and model')

sps = build_sps()
#sps = pread.get_sps(res)
mod = build_model()
thetas = mod.theta_labels()

thetas_50 = []
thetas_16 = []
thetas_84 = []
print('quantiles for all thetas')
for theta in thetas:
    idx = thetas.index(theta)
    chain = [item[idx] for item in res['chain']]
    quan = quantile(chain, [.16, .5, .84])
    thetas_50.append(quan[1])
    thetas_16.append(quan[0])
    thetas_84.append(quan[2])


mod_50 = mod.mean_model(thetas_50, obs, sps)
massfrac_50 = mod_50[-1]
spec = mod_50[0]
mod_16 = mod.mean_model(thetas_16, obs, sps)
massfrac_16 = mod_16[-1]
mod_84 = mod.mean_model(thetas_84, obs, sps)
massfrac_84 = mod_84[-1]

print('mass and Z')
mass = thetas_50[thetas.index('massmet_1')]
mass_50 = thetas_50[thetas.index('massmet_1')]
mass_16 = thetas_16[thetas.index('massmet_1')]
mass_84 = thetas_84[thetas.index('massmet_1')]
Z_50 = thetas_50[thetas.index('massmet_2')]
Z_16 = thetas_16[thetas.index('massmet_2')]
Z_84 = thetas_84[thetas.index('massmet_2')]

print('sfr')
total_mass = 10**mass
time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
print('time: ',len(time_bins_log))
print(time_bins_log)
zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
#print(zfrac_idx)
zfrac_chain = [item[zfrac_idx[0]:zfrac_idx[-1]+1]  for item in res['chain']]
#print(np.shape(zfrac_chain))
sfr_chain = []
for i in range(len(zfrac_chain)):
    sfr_chain.append(zfrac_to_sfr(total_mass, zfrac_chain[i], time_bins_log))
print(np.shape(sfr_chain))
sfr_quan = []
for i in range(np.shape(zfrac_chain)[1]+1): 
    sfr_quan.append(quantile([item[i] for item in sfr_chain], [.16, .5, .84]))

sfr_50.append([item[1] for item in sfr_quan])
sfr_16.append([item[0] for item in sfr_quan])
sfr_84.append([item[2] for item in sfr_quan])
time = []
#for val in time_bins_log:
#    time.append((10**val / 1.0e9)) #Gyr


#print('sfr:',len(sfr_50[0]))
#print(sfr_50[0][::-1])

data = {'Galaxy' : galaxy,  'True Mass' : int_mass, 'Mass_50' : mass_50, 'SFR_50': sfr_50, 'SFR_84': sfr_84, 'SFR_16': sfr_16, 
        'Mass_16' : mass_16, 'Mass_84' : mass_84, 'Z_50' : Z_50, 'Z_16' : Z_16, 
        'Z_84' : Z_84, 'Massfrac' : massfrac_50, 'Spec': spec, 'theta_labels': thetas, 'SPS': sps}


np.savez(prosp_dir+'/galaxy_'+str(galaxy)+'.npz', data=data)
