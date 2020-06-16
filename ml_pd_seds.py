import numpy as np
import pandas as pd
import astropy.constants as constants
import astropy.units as u
from hyperion.model import ModelOutput
import tqdm

#galaxies = pd.read_pickle('/orange/narayanan/s.lower/simba/ml_SEDs_z0.pkl')['ID']

#galaxies = 

pds = glob.glob('/orange/narayanan/s.lower/

spec_list = []
wave_list = []
for galaxy in tqdm.tqdm(galaxies):
    
    m = ModelOutput('/ufrc/narayanan/s.lower//pd_runs/simba_m25n512/snap305/mist_pd/snap305/snap305.galaxy'+"{:03d}".format(galaxy)+".rtout.sed")
    wav,flux = m.get_sed(inclination=0,aperture=-1)
    wav  = np.asarray(wav)*u.micron
    truncate_llim = (np.abs(wav.value - 0.005)).argmin()
    truncate_ulim = (np.abs(wav.value - 1000.)).argmin()
    #print(truncate_llim)
    #print(wav[-1])
    flux = np.asarray(flux)*u.erg/u.s
    dl = (10. * u.pc).to(u.cm)
    flux /= (4.*3.14*dl**2.)
    nu = constants.c.cgs/(wav.to(u.cm))
    nu = nu.to(u.Hz)
    flux /= nu
    flux = flux[truncate_ulim:truncate_llim].to(u.mJy)

    spec_list.append(flux.value)
    wave_list.append(wav[truncate_ulim:truncate_llim].value)

data = {'ID' : galaxies, 'Spectrum [mJy]' : spec_list, 'Wavelength [micron]' : wave_list}

s = pd.DataFrame(data, index=np.arange(len(galaxies)))
s.to_pickle('/orange/narayanan/s.lower/simba/ml_spectra_z0.pkl') 



    
