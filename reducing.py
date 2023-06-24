import h5py
import numpy as np
from util import hdf5lib as h5
import os
import glob
import time

pixels = 256
spx_mass_lim = 10**(-3.0)
spx_sfr_lim = 10**(-10.0)

data_dir = '/n/scratchlfs/vogelsberger/jqi/galaxy/L75n1820TNG/'

idx_list = [99, 91, 84, 78, 72]
#idx_list = [67, 50, 40, 33]

tot_sm, tot_gm, tot_sfr, tot_metal, tot_rad, tot_fg = [], [], [], [], [], []

for snap_idx in idx_list:
	ts = time.time()
	file_list = glob.glob(data_dir + str(snap_idx).zfill(3) + '/rs_maps*')
	#file_list = glob.glob('./test/rs_maps*')
	print file_list
	out_file = data_dir + str(snap_idx).zfill(3) + '/spx_data.hdf5'
	#out_file = './test/spx_data.hdf5'
	for file in file_list:
		print file
		cur_f = h5py.File(file, 'r')
		grps = cur_f.keys()
		for grp_name in grps:
			cur_grp = cur_f[grp_name]
			subs = cur_grp.keys()
			for sub_name in subs:
				cur_sub = cur_grp[sub_name]
				stellar_mass_map = np.array(cur_sub['Stellar_Mass']).flatten()
				gas_SFR_map = np.array(cur_sub['SFR']).flatten()
				gas_metallicity_map = np.array(cur_sub['Gas_Metallicity']).flatten()
				gas_mass_map = np.array(cur_sub['Gas_Mass']).flatten()
				spx_sm, spx_gm, spx_sfr, spx_metal, spx_rad, spx_fg = [], [], [], [], [], []

				for i in xrange(len(stellar_mass_map)):
					if (stellar_mass_map[i] < spx_mass_lim) or \
						(gas_SFR_map[i] < spx_sfr_lim):
						continue
					#rad = (((i / pixels) - pixels/2 + 0.5)**2.0 + ((i % pixels) - pixels / 2 + 0.5)**2.0)**0.5
					
					spx_sm.append(stellar_mass_map[i])
					spx_gm.append(gas_mass_map[i])
					spx_sfr.append(gas_SFR_map[i])
					spx_metal.append(gas_metallicity_map[i])
					#spx_rad.append(rad)
					spx_fg.append(gas_mass_map[i] / (gas_mass_map[i] + stellar_mass_map[i]))

                                        tot_sm.append(stellar_mass_map[i])
                                        tot_gm.append(gas_mass_map[i])
                                        tot_sfr.append(gas_SFR_map[i])
                                        tot_metal.append(gas_metallicity_map[i])
                                        #tot_rad.append(rad)
                                        tot_fg.append(gas_mass_map[i] / (gas_mass_map[i] + stellar_mass_map[i]))


				out = h5py.File(out_file, 'a')
				if not out.__contains__(grp_name):
					out.create_group(grp_name)
				out_grp = out[grp_name]
				if not out_grp.__contains__(sub_name):
					out_grp.create_group(sub_name)
				else:
					continue
				out_sub = out_grp[sub_name]
				out_sub['spx_sm'] = spx_sm
				out_sub['spx_gm'] = spx_gm
				out_sub['spx_sfr'] = spx_sfr
				out_sub['spx_metal'] = spx_metal
				out_sub['spx_fg'] = spx_fg
				#out_sub['spx_rad'] = spx_rad
				#out_sub.attrs = cur_sub.attrs
				out.close()
		cur_f.close()

out = h5py.File(out_file, 'a')
out.create_group('tot')
tot = out['tot']
tot['spx_sm'] = tot_sm
tot['spx_gm'] = tot_gm
tot['spx_sfr'] = tot_sfr
tot['spx_metal'] = tot_metal
tot['spx_fg'] = tot_fg
#tot['spx_rad'] = tot_rad
out.close()
