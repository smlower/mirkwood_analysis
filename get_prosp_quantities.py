from corner import quantile
import prospect.io.read_results as pread
from prospect.models.transforms import zfrac_to_sfrac, zfrac_to_sfr, zfrac_to_masses
import numpy as np
import sys
import glob, os
import pandas as pd
from tqdm.auto import tqdm


#------------------------------------------------- 

def get_best(res, **kwargs):
    imax = np.argmax(res['lnprobability'])
    theta_best = res['chain'][imax, :].copy()
    return theta_best

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx


#-------------------------------------------------



galaxy = sys.argv[1]
sim = sys.argv[2]
snr = sys.argv[3]
z = sys.argv[4]


prosp_dir = '/orange/narayanan/s.lower/prospector/mirkwood_comp/'+str(sim)+'/z'+str(z)+'/snr_'+str(snr)+'/'


print('importing model and sps')
sys.path.append(prosp_dir)
print(prosp_dir)
from run_prosp import build_model, build_sps



sfr_50 = [] 
sfr_16 = []
sfr_84 = []

print('now reading files')

infile = prosp_dir+'/galaxy_'+str(int(galaxy))+'_SNR_'+str(snr)+'.h5'

for prosp_output in glob.glob(infile):
    print(prosp_output)
    res, obs, mod = pread.results_from(prosp_output)

print('building sps and model')

sps = build_sps()
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
spec = mod_50[0]
print('median quantities')
print('    mass and Z')
mass = thetas_50[thetas.index('massmet_1')]
mass_50 = thetas_50[thetas.index('massmet_1')]
mass_16 = thetas_16[thetas.index('massmet_1')]
mass_84 = thetas_84[thetas.index('massmet_1')]
Z_50 = thetas_50[thetas.index('massmet_2')]
Z_16 = thetas_16[thetas.index('massmet_2')]
Z_84 = thetas_84[thetas.index('massmet_2')]

print('    sfr')
total_mass = 10**mass
time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
zfrac_idx = [i for i, s in enumerate(thetas) if 'z_fraction' in s]
zfrac_chain = [item[zfrac_idx[0]:zfrac_idx[-1]+1]  for item in res['chain']]
sfr_chain = []
for i in range(len(zfrac_chain)):
    sfr_chain.append(zfrac_to_sfr(total_mass, zfrac_chain[i], time_bins_log))
sfr_quan = []
for i in range(np.shape(zfrac_chain)[1]+1): 
    sfr_quan.append(quantile([item[i] for item in sfr_chain], [.16, .5, .84]))

sfr_50.append([item[1] for item in sfr_quan])
sfr_16.append([item[0] for item in sfr_quan])
sfr_84.append([item[2] for item in sfr_quan])


print('    dust mass')
print('        maximum liklihood thetas')
theta_max = get_best(res)

total_mass = 10**theta_max[thetas.index('massmet_1')]
time_bins_log = next(item for item in res['model_params'] if item["name"] == "agebins")['init']
z_frac = theta_max[thetas.index('z_fraction_1'):thetas.index('z_fraction_5')+1]
masses = zfrac_to_masses(total_mass, z_frac, time_bins_log)
converted_sfh = sps.convert_sfh(time_bins_log, masses)
model_sp = sps.ssp
model_sp.params['sfh'] = 3
#model_sp.params['zcontinuous'] = 1
model_sp.params['imf_type'] = 2
model_sp.params['zred'] = 0.0
model_sp.set_tabular_sfh(converted_sfh[0], converted_sfh[1])
model_sp.params['add_dust_emission'] = True
model_sp.params['dust_type'] = 5
model_sp.params['dust2'] = theta_max[thetas.index('dust2')]
model_sp.params['dust_index'] = theta_max[thetas.index('dust_index')]
model_sp.params['logzsol'] = theta_max[thetas.index('massmet_2')]
model_sp.params['duste_gamma'] = theta_max[thetas.index('duste_gamma')]
model_sp.params['duste_umin'] = theta_max[thetas.index('duste_umin')]
model_sp.params['duste_qpah'] = 5.86
dust_mass = model_sp.dust_mass




print('done. saving')
data = {'Galaxy' : galaxy, 'mass_50' : mass_50, 'sfr_50': sfr_50[0], 
        'sfr_84': sfr_84[0], 'sfr_16': sfr_16[0],
        'mass_16' : mass_16, 'mass_84' : mass_84, 'logzsol_50' : Z_50, 'logzsol_16' : Z_16, 'dust_mass': dust_mass,'logzsol_84' : Z_84}


np.savez(prosp_dir+'/galaxy_'+str(galaxy)+'_prosp.npz', data=data)
