import h5py
import numpy as np
from util import hdf5lib as h5
import os
import glob
import time

#idx = 0
max_size = 1024
#max_grp = 0

data_dir = '/n/scratchlfs/vogelsberger/jqi/galaxy/L75n1820TNG/'

idx_list = [99, 91, 84, 78, 72, 67, 50, 40, 33]

for snap_idx in idx_list:
	idx = 0
	max_grp = 0

	ts = time.time()
	file_list = glob.glob(data_dir + str(snap_idx).zfill(3) + '/resolved_maps_p*')
	print file_list

	for file in file_list:
		cur_f = h5py.File(file, 'r')
		grp_list = cur_f.keys()
		for grp in grp_list:
			grp_num = int(grp[3:])
			max_grp = max(max_grp, grp_num)
		cur_f.close()
	print max_grp

	for grp in xrange(max_grp):
		cur_grp_name = 'grp' + str(grp)
		for file in file_list:
			cur_file = h5py.File(file, 'r')
			if not (cur_file.__contains__(cur_grp_name)):
				continue
			grp_to_copy = cur_file[cur_grp_name]
			cur_out = data_dir + str(snap_idx).zfill(3) + '/rs_maps.' + str(idx) + '.hdf5'
			if (os.path.exists(cur_out)):
				cur_size = os.stat(cur_out).st_size / float(1024 * 1024)
				if (cur_size >= max_size):
					idx += 1
			new_file_name = data_dir + str(snap_idx).zfill(3) + '/rs_maps.' + str(idx) + '.hdf5'
			new_file = h5py.File(new_file_name, 'a')
			new_file.copy(grp_to_copy, new_file)
			new_file.close()
		cur_file.close()
	te = time.time()
	print te - ts
