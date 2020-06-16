import h5py
import numpy as np
import sys, os
import numpy as np
import glob
import tqdm

snap = int(sys.argv[1])




snap_dir = '/orange/narayanan/s.lower/simba/m100n1024/desika_filtered_snaps/snap0'+str(snap)+'/'
outfile = '/orange/narayanan/s.lower/simba/m100n1024/desika_filtered_snaps/simba_snap'+str(snap)+'_pos.npz'
#snap_dir = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/snap'+str(snap)+'/'
#outfile = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/simba_snap'+str(snap)+'_pos.npz'


pos = {}
ngalaxies = {}

infiles = sorted(glob.glob(snap_dir+'/galaxy_*.hdf5'))

for i in tqdm.tqdm(range(len(infiles))):
    #print('on galaxy ',i,' of ',len(gal))
    try:
        infile = h5py.File(snap_dir+'/galaxy_'+str(i)+'.hdf5', 'r')
    except:
        continue
    pos['galaxy'+str(i)] = {}
    

    gas_masses = infile['PartType0']['Masses']
    gas_coords = infile['PartType0']['Coordinates']
    star_masses = infile['PartType4']['Masses']
    star_coords = infile['PartType4']['Coordinates']
    total_mass = np.sum(gas_masses) + np.sum(star_masses)
    
    x_pos = (np.sum(gas_masses * gas_coords[:,0]) + np.sum(star_masses * star_coords[:,0])) / total_mass 
    y_pos = (np.sum(gas_masses * gas_coords[:,1]) + np.sum(star_masses * star_coords[:,1])) / total_mass
    z_pos = (np.sum(gas_masses * gas_coords[:,2]) + np.sum(star_masses * star_coords[:,2])) / total_mass


    pos['galaxy'+str(i)]['snap'+str(snap)] = np.array([x_pos, y_pos, z_pos])
    infile.close()
#print(x_pos, y_pos, z_pos)
ngalaxies['snap'+str(snap)] = len(infiles)




np.savez(outfile, ngalaxies=ngalaxies, pos=pos)
