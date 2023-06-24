import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
import datetime
from astropy.io import fits

from util import cosmo_tools

import simread.readsubfHDF5 as read_subf
import simread.readhaloHDF5 as rh
import simread.readsnapHDF5 as rs

import visualization.image_maker as im
import h5py
from mpi4py import MPI

from fcs import *
#######################################################
# Input Dataset
data_dir = '/orange/paul.torrey/IllustrisTNG/Runs'
model = '/L75n455TNG/'
snap_idx = 99

out_dir = '/orange/paul.torrey/jqi/Project/TNG_analysis/'
dataset = data_dir + model #+ 'output/'
snapbase = 'snap'

#######################################################
# Property Maps Parameter
galaxy_mass_threshold 	= 10.0**(-1)		# Minimum host galaxy stellar mass to process, in 10**(10) M_sun
galaxy_SFR_threshold 	= 10.0**(-2)		# Minimum host galaxy SFR, in M_sun per year

FOV		= 128.0				# FOV for property maps, in comoving kpc, actual map will be 2FOV * 2FOV
pixels		= 256				# number of pixels in property maps, actual map will be pixels * pixels
						# Resolution of the map will be ( 2 * FOv / pixels )
face_flag	= True				# Rotate galaxies to face-on
h		= 0.6774
non_CGM		= True				# cutting off all CGM gases
center_pixel	= True				# center pixel enclose the halo center

recover		= False				# Recover from a previous stopped run
clip		= False				# Clip galaxy with furthest particles
save		= 'hdf5'			# Output format for property maps, currently supporting 'hdf5' and 'fits'
max_size	= 9999999			# Max size of each output file, in MB

######################################################
# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()
n_tasks = size
this_task = rank
sys.stdout.write(
        "Start work in process %d of %d on %s.\n"
        % (rank, size, name))

######################################################
# Main Logic

# Load in the snapshot
#blockPrint()
HR = rh.HaloReader( dataset, snap_idx, snapbase=snapbase)
cat = read_subf.subfind_catalog(data_dir + model, snap_idx, \
        keysel= ['GroupNsubs', 'SubhaloMass', 'SubhaloSFR', \
        'SubhaloGasMetallicity', 'SubhaloStellarPhotometrics', \
        'SubhaloPos', 'SubhaloLenType', 'SubhaloMassType', 'SubhaloVel'])

#enablePrint()

snapshot = dataset + 'snapdir_' + str(snap_idx).zfill(3) + '/' \
        + snapbase + '_' + str(snap_idx).zfill(3) + '.0.hdf5'
header = rs.snapshot_header(snapshot)
if rank == 0:
        print snapshot
        print 'Loading Snapshot with Redshift ' + str(np.round(1 / header.time - 1, decimals=3))

a = header.time
redshift = 1 / header.time - 1

# Reset FOV to account for h and redshift
FOV = FOV * h * (1 + redshift)

ngroups = HR.cat.ngroups
nsubs = HR.cat.nsubs
nsub_list = cat.GroupNsubs
center_poss = cat.SubhaloPos
center_vels = cat.SubhaloVel

if rank == 0:
        print 'There are ' + str(ngroups) + ' groups, and ' + str(nsubs) + ' subhalos.'

# Start from last failure
halo_idx, grp_start, sub_start = generate_index(recover)

cnt = 0
part = 0
status_idx = 0

for grp_num in xrange(grp_start, len(nsub_list)):
        if (grp_num != grp_start):
                sub_start = 0
        if (grp_num % size) != rank:
                continue
        for sub_num in xrange(sub_start, nsub_list[grp_num]):
                halo_idx = sum(nsub_list[:grp_num]) + sub_num
                print '###################################################################'
                if (check_exist(out_dir, rank, grp_num, sub_num)):
                        print 'Halo ' + str(halo_idx) + ', data already existed'
                        continue

                print 'start reducing ' + str(halo_idx) + ', on process ' + str(rank)
                # Global properties for the subhalo
                halo_mass = cat.SubhaloMassType[halo_idx][4] / h
                halo_sfr = cat.SubhaloSFR[halo_idx]
                halo_metallicity = cat.SubhaloGasMetallicity[halo_idx]
                halo_color = cat.SubhaloStellarPhotometrics[halo_idx]
                print halo_mass, halo_sfr, halo_metallicity

                # Select galaxies with lower limit for mass and SFR
                if (halo_mass < galaxy_mass_threshold) or (halo_sfr < galaxy_SFR_threshold):
                        print 'Halo ' + str(halo_idx) + ', below threshold.'
                        continue

                # Check if there are stellar and gas data
                sp = cat.SubhaloLenType[halo_idx][4]
                gp = cat.SubhaloLenType[halo_idx][0]

                if (sp < 64) or (gp < 64):
                        print 'Halo ' + str(halo_idx) + ', too few particles.'
                        continue
                cnt += 1

                try:
                        blockPrint()

                	stellar_mass_map, gas_SFR_map, gas_metallicity_map, gas_mass_map = \
                                custom_mkim(HR, grp_num, sub_num, center_poss[halo_idx], \
                                        center_vels[halo_idx], FOV=FOV, pixels=pixels, \
					face_flag=face_flag, non_CGM=non_CGM, \
					center_pixel=center_pixel)

                        enablePrint()
                        print 'Maps generated.' + str(rank)
                except:
                        print 'Image_maker failed.'
                        continue

                # Clip the image to reduce storage, just need all pixels with mass higher than threshold
                if clip:
                        x_range, y_range = calculate_range(stellar_mass_map, 10, FOV=FOV, pixels=pixels, h=h)
                        stellar_mass_map = clip_image(stellar_mass_map, y_range=y_range, x_range=x_range)
                        gas_SFR_map = clip_image(gas_SFR_map, y_range=y_range, x_range=x_range)
                        gas_metallicity_map = clip_image(gas_metallicity_map, y_range=y_range, x_range=x_range)
                        gas_mass_map = clip_image(gas_mass_map, y_range=y_range, x_range=x_range)

                hnames = ['grp_num', 'sub_num', 'n_star', 'n_gas', 'mass', 'SFR', 'metal', 'U', 'B', 'V', 'K', 'g', 'r', 'i', 'z']
                hvalues = [grp_num, sub_num, sp, gp, halo_mass, halo_sfr, halo_metallicity, halo_color[0], halo_color[1], halo_color[2], halo_color[3], halo_color[4], halo_color[5], halo_color[6], halo_color[7]]
                image_names = ['Stellar_Mass', 'SFR', 'Gas_Metallicity', 'Gas_Mass']
                images = [stellar_mass_map, gas_SFR_map, gas_metallicity_map, gas_mass_map]

                ######################################3
                if save == 'hdf5':
                        print 'Start saving...'
                        switched = False
                        while not switched:
                                cur_file_name = out_dir + 'resolved_maps_p' + str(rank)+ '.' + str(part) + '.hdf5'
                                part, switched = switch_file(cur_file_name, part, max_size)

                        # Check if current group exists
                        file_list = glob.glob(out_dir + 'resolved_maps_p' + str(rank) + '.*')
                        for file in file_list:
                                fn = file
                                grp_name = 'grp' + str(grp_num)
                                try:
                                        f = h5py.File(fn, 'r')
                                except:
                                        continue
                                if (f.__contains__(grp_name)):
                                        part = int(file.split('.')[-2])
                                f.close()

                        try:
                                save_to_hdf5(out_dir, grp_num, sub_num, part, rank, images=images, image_names=image_names, header_names=hnames, header_values=hvalues)
                        except:
                                print ('Saving Error ########################')
                                continue

