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
caesar_file = glob.glob('/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0'+"{:03d}".format(args.snapnum)+'*.hdf5')
obj = caesar.load(caesar_file[0])

glist = obj.galaxies[int(args.galaxy)].glist
slist = obj.galaxies[int(args.galaxy)].slist


outfile = '/orange/narayanan/s.lower/simba/filtered_snapshots/snap'+"{:03d}".format(args.snapnum)+'/'
#outfile = '/orange/narayanan/s.lower/simba/desika_filtered_snaps/test/'
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


print('done copying attributes, going to edit header now')
outfile_reload = outfile+'galaxy_'+str(args.galaxy)+'.hdf5'

re_out = h5py.File(outfile_reload)                                                  
re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([len(glist), 0, 0, 0, len(slist), 0]))  
re_out['Header'].attrs.modify('NumPart_Total', np.array([len(glist), 0, 0, 0, len(slist), 0]))  

re_out.close()
