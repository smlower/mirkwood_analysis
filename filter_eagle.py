import h5py
import read_eagle
import numpy as np
from mpi4py import MPI
import pandas as pd
import sys, os
import tqdm

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

snap_num = sys.argv[1]
z = sys.argv[2]
galaxy = sys.argv[3]

z_padded = "{:07.3f}".format(float(z))
a = z_padded.split('.')
snap_z = 'z'+a[0]+'p'+a[1]
fname = 'snap_'+"{:03d}".format(int(snap_num))+'_'+snap_z+'.0.hdf5'
dirname = 'snapshot_'+"{:03d}".format(int(snap_num))+'_'+snap_z

IDs = pd.read_csv('/orange/narayanan/s.lower/eagle/filtered_snapshots/galaxy_lists/snap'+snap_num+'_halo_galaxy.csv')
group = IDs['GroupNumber'][int(galaxy)]
subgroup = IDs['SubGroupNumber'][int(galaxy)]

print('Read in galaxy ID')

snap_dir = '/orange/narayanan/s.lower/eagle/m50n752/snapshots/RefL0050N0752/'+dirname
snap = read_eagle.EagleSnapshot(snap_dir+'/'+fname)

snap.select_region(0, 50. * 0.6777, 0, 50. * 0.6777, 0, 50. * 0.6777)
snap.split_selection(comm_rank, comm_size)
print('Selected region')
outdir = '/orange/narayanan/s.lower/eagle/filtered_snapshots/snap023/'
output_file = h5py.File(outdir+'/galaxy_'+str(galaxy)+'.hdf5', 'w')
#output_file = h5py.File('/orange/narayanan/s.lower/eagle/filtered_snapshots/snap'+snap_num+'/galaxy_'+str(galaxy)+'.hdf5', 'w')
input1 = h5py.File(snap_dir+'/'+fname, 'r')
#output_file.copy(input1['Header'], 'Header')

attrs = ['Config',
 'Constants',
 'HashTable',
 'Header',
 'Parameters',
 'RuntimePars',
 'Units']

for key in attrs:

    output_file.copy(input1[key], key)

input1.close()

print('Copied Header and others')

gas_attr = snap.datasets(0)
print('Found gas attributes')
star_attr = snap.datasets(4)
print('Found star attributes')
gas_groups = snap.read_dataset(0, 'GroupNumber')
print('got gas groups')
gas_subgroups = snap.read_dataset(0, 'SubGroupNumber')
print('got gas subgroups')
star_groups = snap.read_dataset(4, 'GroupNumber')
print('got star groups')
star_subgroups = snap.read_dataset(4, 'SubGroupNumber')
print('got star subgroups')

mask_g = np.logical_and(gas_groups == group, gas_subgroups == subgroup)
print('got gas mask')
mask_s = np.logical_and(star_groups == group, star_subgroups == subgroup)
print('got star mask')
output_file.create_group('PartType0')
output_file.create_group('PartType4')

print('going into gas attributes')
for attr in tqdm.tqdm(gas_attr):
    output_file['PartType0'][attr[1:]] = snap.read_dataset(0, attr)[mask_g]
    if attr[1:] == 'Mass':
        output_file['PartType0']['Masses'] = snap.read_dataset(0, attr)[mask_g]
print('going into star attributes')
for attr in tqdm.tqdm(star_attr):
    output_file['PartType4'][attr[1:]] = snap.read_dataset(4, attr)[mask_s]
    if attr[1:] == 'Mass':
        output_file['PartType4']['Masses'] = snap.read_dataset(4, attr)[mask_s]
output_file.close()




re_out = h5py.File(outdir+'/galaxy_'+str(galaxy)+'.hdf5', 'r+')

slist = len(re_out['PartType4']['Coordinates'])
glist = len(re_out['PartType0']['Coordinates'])
re_out['Header'].attrs.modify('NumPart_ThisFile', np.array([glist, 0, 0, 0, slist, 0]))
re_out['Header'].attrs.modify('NumPart_Total', np.array([glist, 0, 0, 0, slist, 0]))
re_out['Header'].attrs.modify('NumFilesPerSnapshot', 1)
re_out.close()

    

