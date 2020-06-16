import h5py
import caesar
import argparse, sys
import glob
import numpy as np
import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('snapnum', type=int)
parser.add_argument('galaxy', type=int)
args = parser.parse_args()

ds = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_'+"{:03d}".format(args.snapnum)+'.hdf5'
#ds = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m100n1024/snap_m100n1024_078.hdf5'
#caesar_file = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m100n1024/caesar_m100n1024_078.hdf5'
caesar_file = glob.glob('/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0'+"{:03d}".format(args.snapnum)+'*.hdf5')
obj = caesar.quick_load(caesar_file[0])

if obj.galaxies[int(args.galaxy)].masses['gas'] < 1.:
    sys.exit()

ids = np.load('/orange/narayanan/s.lower/simba/snap'+str(args.snapnum)+'_galaxy.npz', allow_pickle=True)
gal = ids['galid'][args.galaxy]

glist = obj.galaxies[int(args.galaxy)].glist
slist = obj.galaxies[int(args.galaxy)].slist
#bhlist = obj.galaxies[int(args.galaxy)].bhlist
#parent_halo = obj.galaxies[int(args.galaxy)].parent_halo_index
#dmlist = obj.halos[parent_halo].dmlist


outfile = '/orange/narayanan/s.lower/simba/filtered_snapshots/snap'+str(args.snapnum)+'/'
#outfile = '/orange/narayanan/s.lower/simba/m100n1024/desika_filtered_snaps/snap078/'
with h5py.File(ds, 'r') as input_file, h5py.File(outfile+'galaxy_'+str(args.galaxy)+'.hdf5', 'w') as output_file:
    output_file.copy(input_file['Header'], 'Header')
    print('starting with gas attributes now')
    output_file.create_group('PartType0')
    for k in tqdm.tqdm(input_file['PartType0']):
        output_file['PartType0'][k] = input_file['PartType0'][k][:][glist]
    print('moving to star attributes now')
    output_file.create_group('PartType4')
    for k in tqdm.tqdm(input_file['PartType4']):
        output_file['PartType4'][k] = input_file['PartType4'][k][:][slist]

    #output_file.create_group('PartType5')
    #for k in input_file['PartType5']:
    #    output_file['PartType5'][k] = input_file['PartType5'][k][:][bhlist]
        

    #output_file.create_group('PartType1')
    #for k in input_file['PartType1']:
    #    output_file['PartType1'][k] = input_file['PartType1'][k][:][dmlist]


print('done copying attributes, going to edit header now')
outfile_reload = outfile+'galaxy_'+str(args.galaxy)+'.hdf5'

re_out = h5py.File(outfile_reload)                                                  
re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([len(glist), 0, 0, 0, len(slist), 0]))  
re_out['Header'].attrs.modify('NumPart_Total', np.array([len(glist), 0, 0, 0, len(slist), 0]))  

re_out.close()
