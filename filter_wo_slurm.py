import h5py
import caesar
import argparse, sys
import glob
import numpy as np
from tqdm.auto import tqdm
import multiprocessing


def filter_galaxy(galaxy):
    glist = obj.galaxies[int(galaxy)].glist
    slist = obj.galaxies[int(galaxy)].slist
    outfile = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/test/'
    #outfile = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/snap'+"{:03d}".format(snap)+'/'
    with h5py.File(ds, 'r') as input_file, h5py.File(outfile+'galaxy_'+str(galaxy)+'.hdf5', 'w') as output_file:
        output_file.copy(input_file['Header'], 'Header')
        #print('starting with gas attributes now')
        output_file.create_group('PartType0')
        for k in input_file['PartType0']:
            output_file['PartType0'][k] = input_file['PartType0'][k][:][glist]
        #print('moving to star attributes now')
        output_file.create_group('PartType4')
        for k in input_file['PartType4']:
            output_file['PartType4'][k] = input_file['PartType4'][k][:][slist]                                                                                                    

    #print('done with galaxy',galaxy)
    #print('done copying attributes, going to edit header now')
    outfile_reload = outfile+'galaxy_'+str(galaxy)+'.hdf5'
        
    re_out = h5py.File(outfile_reload)
    re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([len(glist), 0, 0, 0, len(slist), 0]))
    re_out['Header'].attrs.modify('NumPart_Total', np.array([len(glist), 0, 0, 0, len(slist), 0]))
        
    re_out.close()
    print('done with galaxy',galaxy)
    
    return 



snap = 183
ds = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_'+"{:03d}".format(snap)+'.hdf5'
caesar_file = glob.glob('/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0'+"{:03d}".format(snap)+'*.hdf5')
obj = caesar.load(caesar_file[0])


if __name__ == '__main__':
    with multiprocessing.Pool(16) as pool:
        pool.map(filter_galaxy, range(0, 2000))
        pool.close()
        pool.join()
