import h5py
import numpy as np
import sys, os
import tqdm


galaxy = int(sys.argv[1])

base_path = '/orange/narayanan/s.lower/TNG/'
snap_num = 95
outfile = base_path+'/filtered_snapshots/snap095/galaxy_'+str(galaxy)+'.hdf5'


gal_list = np.load('/orange/narayanan/s.lower/TNG/galaxy_selections/tng_galaxy_selection_z05.npz')['ID']
gal = gal_list[int(galaxy)]

snap1 = base_path+'/output/snapdir_095/snap_095.0.hdf5'
fullfile = h5py.File(snap1, 'r')

print('creating snapshot for galaxy '+str(gal))

with h5py.File(base_path+'simulation.hdf5','r') as input_file, h5py.File(outfile, 'w') as output_file: 
    output_file.copy(fullfile['Header'], 'Header')
    output_file.copy(fullfile['Config'], 'Config')
    
    print('copying gas attributes')
    output_file.create_group('PartType0')
    start = input_file['/Offsets/'+str(snap_num)+'/Subhalo/SnapByType'][gal, 0]
    glength = input_file['/Groups/'+str(snap_num)+'/Subhalo/SubhaloLenType'][gal, 0]
    print('gas start: ',start)
    print('gas end: ',start+glength)
    for k in tqdm.tqdm(input_file['/Snapshots/'+str(snap_num)+'/PartType0/']):
        output_file['PartType0'][str(k)] = input_file['/Snapshots/'+str(snap_num)+'/PartType0/'+str(k)][start:start+glength]
    
    print('copying star attributes')
    output_file.create_group('PartType4')
    start = input_file['/Offsets/'+str(snap_num)+'/Subhalo/SnapByType'][gal, 4]
    slength = input_file['/Groups/'+str(snap_num)+'/Subhalo/SubhaloLenType'][gal, 4]
    print('star start: ',start)
    print('star end: ',start+slength)
    for k in tqdm.tqdm(input_file['/Snapshots/'+str(snap_num)+'/PartType4/']):
        c = str(k)
        if c == 'StellarAssembly':
            continue
        output_file['PartType4'][str(k)] = input_file['/Snapshots/'+str(snap_num)+'/PartType4/'+str(k)][start:start+slength]


re_out = h5py.File(outfile, 'r+')

re_out['Header'].attrs.modify('NumFilesPerSnapshot', 1)
re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([glength, 0, 0, 0, slength, 0]))
re_out['Header'].attrs.modify('NumPart_Total', np.array([glength, 0, 0, 0, slength, 0]))

re_out.close()

print('galaxy '+str(galaxy)+' is done.')
