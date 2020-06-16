import yt
import numpy as np
import pandas as pd
import fsps
from hyperion.model import ModelOutput
from sedpy.observate import load_filters
import glob, tqdm, sys
import astropy.constants as constants
import astropy.units as u


snap_num = str(sys.argv[1])
z = str(sys.argv[2])
def ssfr_relation(mstar):
    return (10**-13) * (10**mstar)

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx


galex = ['galex_FUV', 'galex_NUV']
hst_wfc3_uv  = ['wfc3_uvis_f275w', 'wfc3_uvis_f336w', 'wfc3_uvis_f475w','wfc3_uvis_f555w', 'wfc3_uvis_f606w', 'wfc3_uvis_f814w']
sdss = ['sdss_i0']
hst_wfc3_ir = ['wfc3_ir_f105w', 'wfc3_ir_f125w', 'wfc3_ir_f140w', 'wfc3_ir_f160w']
spitzer_irac = ['spitzer_irac_ch1']
spitzer_mips = ['spitzer_mips_24']
herschel_pacs = ['herschel_pacs_70', 'herschel_pacs_100', 'herschel_pacs_160']
herschel_spire = ['herschel_spire_250', 'herschel_spire_350', 'herschel_spire_500']
jwst_miri = ['jwst_f560w', 'jwst_f770w', 'jwst_f1000w', 'jwst_f1130w', 'jwst_f1280w', 'jwst_f1500w', 'jwst_f1800w']
jwst_nircam = ['jwst_f070w', 'jwst_f090w', 'jwst_f115w', 'jwst_f150w', 'jwst_f200w', 'jwst_f277w']
filternames = (galex + hst_wfc3_uv +  hst_wfc3_ir + jwst_miri + jwst_nircam + spitzer_irac + spitzer_mips + herschel_pacs + herschel_spire)
filters_unsorted = load_filters(filternames)
waves_unsorted = [x.wave_mean for x in filters_unsorted]
filters = [x for _,x in sorted(zip(waves_unsorted,filters_unsorted))]

#print('Loading SP')
#fsps_ssp = fsps.StellarPopulation(sfh=0,
#                zcontinuous=1,
#                imf_type=2,
#                zred=0.0, add_dust_emission=False)
#solar_Z = 0.0196


snaps_dir = '/orange/narayanan/s.lower/simba/filtered_snapshots/snap'+snap_num+'/'
pd_dir = '/orange/narayanan/s.lower/simba/pd_runs/snap'+snap_num+'/'
snaps = glob.glob(snaps_dir+'/galaxy*.hdf5')
pds = glob.glob(pd_dir+'/*.rtout.sed')

print('initializing lists')
flux_Jy = []
fluxe_Jy = []
stellar_mass_Msun = []
dust_mass_Msun = []
sfr100 = []
metallicity_logzsol = []
gal_count = []
filter_list = []

for galaxy in tqdm.tqdm(range(len(pds))):
    
    try:
        m = ModelOutput(pds[galaxy])
        #ds = yt.load(snaps_dir+'/galaxy_'+str(galaxy)+'.hdf5')
        wav, flux = m.get_sed(inclination=0, aperture=-1)
   
    except:
        continue
    
    gal = pds[galaxy].split('galaxy')[1].split('.')[0]
    gal_count.append(int(gal))
    
    #get mock photometry
    #wav, flux = m.get_sed(inclination=0, aperture=-1)
    wav  = np.asarray(wav)*u.micron #wav is in micron                                                                                                        
    wav = wav.to(u.AA)
    flux = np.asarray(flux)*u.erg/u.s
    dl = (10. * u.pc).to(u.cm)
    flux /= (4.*3.14*dl**2.)
    nu = constants.c.cgs/(wav.to(u.cm))
    nu = nu.to(u.Hz)
    flux /= nu
    flux = flux.to(u.Jy)

    flx = []
    flxe = []
    for i in range(len(filters)):
        flux_range = []
        wav_range = []
        for j in filters[i].wavelength:
            flux_range.append(flux[find_nearest(wav.value,j)].value)
            wav_range.append(wav[find_nearest(wav.value,j)].value)
        a = np.trapz(wav_range * filters[i].transmission* flux_range, wav_range, axis=-1)
        b = np.trapz(wav_range * filters[i].transmission, wav_range)
        flx.append(a/b)
        flxe.append(0.03* flx[i])

    flux_Jy.append(np.asarray(flx))
    fluxe_Jy.append(np.asarray(flxe))


    #get properties: M*, SFR, Z, Mdust
    '''print('getting properties')
    ad = ds.all_data()
    star_masses = ad[('PartType4', 'Masses')].in_units('Msun')
    stellar_mass_Msun.append([np.sum(star_masses.value)])
    print('got stellar mass')
    scalefactor = ad[("PartType4", "StellarFormationTime")]
    star_metal = ad[("PartType4", 'metallicity')]
    metallicity_logzsol.append([np.log10(np.sum(star_metal) / solar_Z)])
    print('got stellar metallicity')
    formation_z = (1.0 / scalefactor) - 1.0
    stellar_formation_age = ds.cosmology.t_from_z(formation_z).in_units("Gyr")
    simtime = ds.cosmology.t_from_z(ds.current_redshift).in_units("Gyr")
    stellar_ages = (simtime - stellar_formation_age).in_units("Gyr")

    w50 = np.where(stellar_ages.in_units('Myr').value < 100)[0]
    if len(w50) == 0:
        print('no stars born within the last 100 Myr. Setting sfr according to sSFR relation')
        sfr100.append(ssfr_relation(np.log10(stellar_mass_Msun)))

    else:
        initial_mass = 0.0
        print('getting initial SP mass')
        for star in tqdm.tqdm(w50):
            current_mass = star_masses[star]
            #print('got current mass')
            fsps_ssp.params["logzsol"] = np.log10(star_metal[star] / solar_Z)
            #print('got metallicity')
            mass_remaining = fsps_ssp.stellar_mass
            #print('got surviving mass')
            #print('finding initial mass')
            initial_mass += current_mass / np.interp(np.log10(stellar_ages[star]*1.e9),fsps_ssp.ssp_ages,mass_remaining)  
        
        sfr_50myr = initial_mass/100.e6
        sfr100.append([sfr_50myr])
    
    print('got SFR')
    
    #metal_mass = (ad[("PartType0", "Masses")])*(ad["PartType0", "Metallicity"].value)
    #dust_masses = metal_mass * 0.4
    dust_masses = ad.ds.arr(ad[("PartType0", "Dust_Masses")].value, 'code_mass')
    #dust_masses = ad[("PartType0", "Dust_Masses")]
    dust_mass_Msun.append([np.sum(dust_masses.in_units('Msun').value)])
    print('got dust mass')'''
    filter_list.append([x.name for x in filters])


print('saving SEDs')
sed_data = {'Filters' : filter_list, 'Flux [Jy]' : flux_Jy, 'Flux Err': fluxe_Jy}

s = pd.DataFrame(sed_data, index=np.arange(len(gal_count)))
s.to_pickle('/orange/narayanan/s.lower/simba/simba_ml_SEDs_z'+z+'_blind.pkl')

np.savez('/orange/narayanan/s.lower/simba/blind_count_order_snap305.npz', index=np.arange(len(gal_count)), galaxy_ID=gal_count)

#print('index: ',np.shape(np.arange(len(gal_count))))
#print('flux: ',np.shape(flux_Jy))
#print('filters: ',np.shape([x.name for x in filters]))
#print('galaxies: ',np.shape(gal_count))
#print('unc: ',np.shape(fluxe_Jy))

#print('saving properties')
#prop_data = {'ID' : gal_count, 'stellar_mass' : stellar_mass_Msun, 'dust_mass' : dust_mass_Msun, 'sfr' : sfr100, 'metallicity' : metallicity_logzsol}
#t = pd.DataFrame(prop_data, index = np.arange(len(gal_count)))
#t.to_pickle('/orange/narayanan/s.lower/simba/simba_ml_props_z'+z+'.pkl')


