import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
import datetime
from util import cosmo_tools
import util.calc_hsml as calc_hsml

import simread.readsubfHDF5 as read_subf
import simread.readhaloHDF5 as rh
import simread.readsnapHDF5 as rs

import visualization.image_maker as im
import visualization.contour_makepic as cmakepic
import h5py
from mpi4py import MPI
from util.calc_hsml import get_gas_density_around_stars
import units.springel_units as units

##################################################################################
def center_halo(gas_pos, gas_vel, star_pos, center_pos, center_vel, box_size):
        sx, sy, sz = star_pos[:, 0], star_pos[:, 1], star_pos[:, 2]
        gx, gy, gz = gas_pos[:, 0], gas_pos[:, 1], gas_pos[:, 2]
	vx, vy, vz = gas_vel[:, 0], gas_vel[:, 1], gas_vel[:, 2]

        gx -= center_pos[0]
        gy -= center_pos[1]
        gz -= center_pos[2]
        sx -= center_pos[0]
        sy -= center_pos[1]
        sz -= center_pos[2]
        vx -= center_vel[0]
        vy -= center_vel[1]
        vz -= center_vel[2]

        for i in xrange(len(sx)):
                if sx[i] > box_size / 2: sx[i] -= box_size
                elif sx[i] < -box_size / 2: sx[i] += box_size
                if sy[i] > box_size / 2: sy[i] -= box_size
                elif sy[i] < -box_size / 2: sy[i] += box_size
                if sz[i] > box_size / 2: sz[i] -= box_size
                elif sz[i] < -box_size / 2: sz[i] += box_size

        for i in xrange(len(gx)):
                if gx[i] > box_size / 2: gx[i] -= box_size
                elif gx[i] < -box_size / 2: gx[i] += box_size
                if gy[i] > box_size / 2: gy[i] -= box_size
                elif gy[i] < -box_size / 2: gy[i] += box_size
                if gz[i] > box_size / 2: gz[i] -= box_size
                elif gz[i] < -box_size / 2: gz[i] += box_size

	return gx, gy, gz, vx, vy, vz, sx, sy, sz

def calc_angles(gx, gy, gz, vx, vy, vz, gas_mass, gas_sfr, r_within=5.0):
        r = np.sqrt( gx**2.0 + gy**(2.0) + gz**(2.0) )
        sel_idx = ( r < r_within )
        lz = np.sum( gas_mass * (gx * vy - gy * vx ) )
        lx = np.sum( gas_mass * (gy * vz - gz * vy ) )
        ly = np.sum( gas_mass * (gz * vx - gx * vz ) )
        phi = np.arctan2( ly, lx )
        theta = np.arctan2( np.sqrt(lx**2 + ly**2), lz )

        return phi, theta

def rotate_halo(phi, theta, gx, gy, gz, vx, vy, vz, sx, sy, sz):
        sx_ = -sz * np.sin(theta) + (sx * np.cos(phi) + sy * np.sin(phi) ) * np.cos(theta)
        sy_ = -sx * np.sin(phi)   + sy * np.cos(phi)
        sz_ =  sz * np.cos(theta) + (sx * np.cos(phi) + sy * np.sin(phi) ) * np.sin(theta)

        gx_ = -gz * np.sin(theta) + (gx * np.cos(phi) + gy * np.sin(phi) ) * np.cos(theta)
        gy_ = -gx * np.sin(phi)   + gy * np.cos(phi)
        gz_ =  gz * np.cos(theta) + (gx * np.cos(phi) + gy * np.sin(phi) ) * np.sin(theta)

        vx_ = -vz * np.sin(theta) + (vx * np.cos(phi) + vy * np.sin(phi) ) * np.cos(theta)
        vy_ = -vx * np.sin(phi)   + vy * np.cos(phi)
        vz_ =  vz * np.cos(theta) + (vx * np.cos(phi) + vy * np.sin(phi) ) * np.sin(theta)

        return gx_, gy_, gz_, vx_, vy_, vz_, sx_, sy_, sz_

def custom_mkim(HR, grp_num, sub_num, cpos, cvel, box_size=75000.0, FOV=128, pixels=256, little_h=0.6774, non_CGM=True, face_flag=True, center_pixel=True):
        # load data from halo_reader
        star_pos 	= HR.read('POS ', 4, grp_num, sub_num)
        star_mass 	= HR.read('MASS', 4, grp_num, sub_num)
        gas_vel 	= HR.read('VEL ', 0, grp_num, sub_num)
        gas_pos 	= HR.read('POS ', 0, grp_num, sub_num)
        gas_mass 	= HR.read('MASS', 0, grp_num, sub_num)
        gas_sfr 	= HR.read('SFR ', 0, grp_num, sub_num)
        gas_metal 	= HR.read('GZ  ', 0, grp_num, sub_num)
        gas_rho 	= HR.read('RHO ', 0, grp_num, sub_num)
        gas_u           = HR.read('U   ', 0, grp_num, sub_num)
        gas_ne          = HR.read('NE  ', 0, grp_num, sub_num)

	# Cutting off Non-CGM gases
	if non_CGM:
		gas_temp = units.gas_code_to_temperature( gas_u, gas_ne )
		non_cgm = np.log10(gas_temp) < (6.0 + 0.25 * np.log10(gas_rho))
		gas_vel 	= gas_vel[non_cgm]
		gas_pos		= gas_pos[non_cgm]
		gas_mass	= gas_mass[non_cgm]
		gas_sfr		= gas_sfr[non_cgm]
		gas_metal	= gas_metal[non_cgm]
		gas_rho		= gas_rho[non_cgm]

        gh = gas_mass / gas_rho
        #gh = 2.0 * 4.0 * gh ** 0.333 / 5.0
	gh = ( 3.0 * gas_mass / (4.0 * np.pi * gas_rho )  )**0.3333

	# Centering
	gx, gy, gz, vx, vy, vz, sx, sy, sz = \
		center_halo(gas_pos, gas_vel, star_pos, cpos, cvel, box_size)
	# Rotating
	if face_flag:
		phi, theta = calc_angles(gx, gy, gz, vx, vy, vz, gas_mass, gas_sfr)
		gx, gy, gz, vx, vy, vz, sx, sy, sz = \
			rotate_halo(phi, theta, gx, gy, gz, vx, vy, vz, sx, sy, sz)

	# Select star-forming gas
	sf_idx = gas_sfr > 0.0
	gx 		= gx[sf_idx]
	gy 		= gy[sf_idx]
	gz 		= gz[sf_idx]
	gas_mass 	= gas_mass[sf_idx]
	gas_sfr 	= gas_sfr[sf_idx]
	gas_metal 	= gas_metal[sf_idx]
	gh		= gh[sf_idx]

	gx, gy, gz = gx / little_h, gy / little_h, gz / little_h
	sx, sy, sz = sx / little_h, sy / little_h, sz / little_h

	if center_pixel:
		offset = FOV / pixels
		gx += offset
		gy += offset
		sx += offset
		sy += offset


	# Create Maps
	sh = calc_hsml.get_particle_hsml(sx, sy, sz)
	print 'Start Mapping'

        ### stellar-mass
        star_mass_map, image = \
                cmakepic.simple_makepic(sx, sy, weights=star_mass, \
                        hsml=sh, xrange=[-FOV, FOV], yrange=[-FOV, FOV], \
                        pixels=pixels)

        ### gas-mass
        gas_mass_map, image = \
                cmakepic.simple_makepic(gx, gy, weights=gas_mass, \
                        hsml=gh, xrange=[-FOV, FOV], yrange=[-FOV, FOV], \
                        pixels=pixels)

        ### gas-sfr
        gas_sfr_map, image = \
                cmakepic.simple_makepic(gx, gy, weights=gas_sfr, \
                        hsml=gh, xrange=[-FOV, FOV], yrange=[-FOV, FOV], \
                        pixels=pixels)

        ### gas-metallicity
        gas_z_map, image = \
                cmakepic.simple_makepic(gx, gy, weights=(gas_metal * gas_mass), \
                        hsml=gh, xrange=[-FOV, FOV], yrange=[-FOV, FOV], \
                        pixels=pixels)

        gas_z_map = gas_z_map / gas_mass_map / 0.0134
        gas_sfr_map = gas_sfr_map
        star_mass_map = star_mass_map * 1e10 / 1e6 / little_h
        gas_mass_map = gas_mass_map * 1e10 / 1e6 / little_h

        return star_mass_map, gas_sfr_map, gas_z_map, gas_mass_map
#######################################################################################

#######################################################
# Disable Print
def blockPrint():
        sys.stdout = open(os.devnull, 'w')

# Restore Print
def enablePrint():
        sys.stdout = sys.__stdout__

# Read halo_idx, grp_num, sub_num from the checkpoint of last failure
def generate_index(recover):
        # order is halo_idx, grp_num, sub_num
        if recover:
                idx_list = np.loadtxt('log.txt')
                return int(idx_list[0]), int(idx_list[1]), int(idx_list[2])
        else:
                return 0, 0, 0

def save_log(filename='log1.txt'):
        log = open(filename, 'w')
        log.write(str(halo_idx))
        log.write('\n')
        log.write(str(grp_num))
        log.write('\n')
        log.write(str(sub_num))
        log.close()
#######################################################
#######################################################
# Clip resolved maps to desired surface density limit
def calculate_range(mass_map, threshold, FOV=128, pixels=256, h=0.6774):
        reso = 2 * FOV / (pixels * h)
        spx_limit = threshold / reso**(2.0)

        x_low = float('inf')
        x_high = float('-inf')
        y_low = float('inf')
        y_high = float('-inf')

        for i in xrange(len(mass_map)):
                for j in xrange(len(mass_map[i])):
                        if mass_map[i][j] > spx_limit:
                                x_low = min(j, x_low)
                                x_high = max(j, x_high)
                                y_low = min(i, y_low)
                                y_high = max(i, y_high)
        return [x_low, x_high], [y_low, y_high]

def clip_image(image, y_range=[0, 0], x_range=[0, 0]):
        clipped = []
        for y in xrange(y_range[0], y_range[1] + 1):
                clipped.append([])
                for x in xrange(x_range[0], x_range[1] + 1):
                        clipped[len(clipped) - 1].append(image[y][x])
        return clipped
#######################################################
#######################################################
# Save Resolved Maps to Files (fits, hdf5)
def save_to_fits(file_name,  images=[], image_names=[], header_names=[], header_values=[]):
        hdr = fits.Header()
        default_names = ['Steallar_Mass', 'SFR', 'Gas_Metallicity', 'Gas_Mass']
        for hname, hvalue in zip(header_names, header_values):
                hdr[hname] = hvalue
        front = fits.PrimaryHDU([[]], header = hdr)

        if len(image_names) == 0:
                image_names = default_names
        hdu_list = [front]
        for iname, image in zip(image_names, images):
                hdu_list.append(fits.ImageHDU(image, name=iname))

        hdul = fits.HDUList(hdu_list)
        hdul.writeto(file_name, overwrite=True)
        print ('Images Saved to ' + file_name)

def save_to_hdf5(out_dir, grp_num, sub_num, part, rank_num, images=[], image_names=[], header_names=[], header_values=[]):
        file_name = out_dir + 'resolved_maps_p' + str(rank_num) + '.' + str(part) + '.hdf5'
        f = h5py.File(file_name, 'a')
        grp_name = 'grp' + str(grp_num)
        sub_name = 'sub' + str(sub_num)
        if not (f.__contains__(grp_name)):
                f.create_group(grp_name)
        cur_grp = f[grp_name]

        if not (cur_grp.__contains__(sub_name)):
                cur_grp.create_group(sub_name)
        cur_sub = cur_grp[sub_name]

        for hname, hvalue in zip(header_names, header_values):
                cur_sub.attrs[hname] = hvalue

        for iname, image in zip(image_names, images):
                cur_sub[iname] = image

        #print ('Group ' + str(grp_num) + ', Subhalo ' + str(sub_num) + ' : ' + 'Images Saved to ' + file_name)
        f.close()

def switch_file(cur_file_name, part, max_size):
        # max_size in MB
        #cur_file_name = out_dir + 'resolved_maps.' + str(part) + '.hdf5'
        if os.path.exists(cur_file_name):
                cur_size = os.stat(cur_file_name).st_size / float(1024 * 1024)
                if (cur_size >= max_size):
                        return part + 1, False
        return part, True
#######################################################
#######################################################
def check_exist(h5dir, rank_num, grp_num, sub_num):
        file_list = glob.glob(h5dir + 'resolved_maps_p' + str(rank_num) + '.*')
        grp_name = 'grp' + str(grp_num)
        sub_name = 'sub' + str(sub_num)
        for file in file_list:
                file_name = file
                try:
                        f = h5py.File(file_name, 'r')
                except:
                        continue
                if (f.__contains__(grp_name)):
                        g = f[grp_name]
                        if (g.__contains__(sub_name)):
                                f.close()
                                return True
                f.close()
        return False

#######################################################
def write_status(status_dir, idx, ts_start, ts_end, grp_num, cur_subs, nsubs, halos):
        if (grp_num < 2000) or (grp_num <20000 and grp_num % 10 == 0) or (grp_num % 100):
                status_file = status_dir + 'status' + str(idx).zfill(3) + '.txt'
                if grp_num == 0:
                        f = open(status_file, 'w')
                else:
                        f = open(status_file, 'a')

                st_start = datetime.datetime.fromtimestamp(ts_start).strftime('%Y-%m-%d %H:%M:%S')
                st_end = datetime.datetime.fromtimestamp(ts_end).strftime('%Y-%m-%d %H:%M:%S')
                f.write('Finished Reducing Group ' + str(grp_num) + '\n')
                f.write('Number of Subhalos: ' + str(cur_subs) + '\n')
                f.write('Start Time: ' + st_start + '\n')
                f.write('End Time: ' + st_end + '\n')
                f.write('Finished: ' + str(halos/float(nsubs)) + ' of total halos' + '\n')
                f.write('----------------------------------------' + '\n')
                f.close()

