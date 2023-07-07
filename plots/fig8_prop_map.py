import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py
from matplotlib.colors import LogNorm
from matplotlib import ticker
#######################################################
def clip_image(arr, fov, center):
        if len(arr) == 0: return arr
        x_max = min(center[0] + fov, len(arr))
        x_min = max(center[0] - fov, 0)
        y_max = min(center[1] + fov, len(arr[0]))
        y_min = max(center[1] - fov, 0)
	#print 2*fov, x_max, x_min, y_max, y_min
        new_im = np.zeros((2*fov, 2*fov))
        for i in xrange(y_min, y_max):
                for j in xrange(x_min, x_max):
		        new_im[i - y_min][j - x_min] = arr[i][j]
        return new_im

#######################################################
data_dir = '../fixed/'

z = 0.0
FOV = 128
pixels = 256
face_flag = True
h = 0.6774

box_size = 75000.0

n_gal = -1

file = h5py.File('../fixed/rs_maps.0.hdf5', 'r')

sub = file['grp13']['sub0']

print sub.keys()

sm = np.array(sub['Stellar_Mass'])
gm = np.array(sub['Gas_Mass'])
sfr = np.array(sub['SFR'])
z = np.array(sub['Gas_Metallicity'])

sm  = clip_image(sm,  50, [128, 128])
gm  = clip_image(gm,  50, [128, 128])
sfr = clip_image(sfr, 50, [128, 128])
z   = clip_image(z,   50, [128, 128])

sm  = np.log10(sm * 10**(6.0) / h)
gm  = np.log10(gm * 10**(6.0) / h)
sfr = np.log10(sfr / h)
z   = np.log10(z)
############################################
fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.10, 0.10, 0.85, 0.85])
plt.imshow(sm, vmin=7.0, vmax=10.0, cmap='coolwarm')

cbar = plt.colorbar(pad=0.025, orientation='horizontal', fraction=0.08, ticks=[7, 8, 9, 10])
cbar.set_label('$\mathrm{log}(\mathrm{\Sigma}_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2}) $', fontsize=60, labelpad=1.1)
cbar.ax.tick_params(labelsize=60)

labels = ax.get_yticklabels()
labels[0] = ''
ax.set_yticklabels(labels)
labels = ax.get_xticklabels()
labels[0] = ''
ax.set_xticklabels(labels)

plt.savefig('sm_map.pdf')

###########################################
fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.10, 0.10, 0.85, 0.85])
plt.imshow(gm, vmin=7.0, vmax=10.0, cmap='coolwarm')

cbar = plt.colorbar(pad=0.025, orientation='horizontal', fraction=0.08, ticks=[7, 8, 9, 10])
cbar.set_label('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{gas}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2}) $', fontsize=60, labelpad=1.1)
cbar.ax.tick_params(labelsize=60)

labels = ax.get_yticklabels()
labels[0] = ''
ax.set_yticklabels(labels)
labels = ax.get_xticklabels()
labels[0] = ''
ax.set_xticklabels(labels)

plt.savefig('gm_map.pdf')

############################################
fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.10, 0.10, 0.85, 0.85])
plt.imshow(sfr, vmin=-4.0, vmax=0.0)

cbar = plt.colorbar(pad=0.025, orientation='horizontal', fraction=0.08, ticks=[0, -1, -2, -3, -4])
cbar.set_label('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2}\, \mathrm{yr}^{-1}) $', fontsize=60, labelpad=1.1)
cbar.ax.tick_params(labelsize=60)

labels = ax.get_yticklabels()
labels[0] = ''
ax.set_yticklabels(labels)
labels = ax.get_xticklabels()
labels[0] = ''
ax.set_xticklabels(labels)

plt.savefig('sfr_map.pdf')

############################################
fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.10, 0.10, 0.85, 0.85])
plt.imshow(z, vmin=-1, vmax=0.1, cmap='gist_heat')

cbar = plt.colorbar(pad=0.025, orientation='horizontal', fraction=0.08, ticks=[0.0, -0.5, -1.0])
cbar.set_label('$\mathrm{log}(\mathrm{Z} / \mathrm{Z}_{\odot}) $', fontsize=60, labelpad=1.1)
cbar.ax.tick_params(labelsize=60)

labels = ax.get_yticklabels()
labels[0] = ''
ax.set_yticklabels(labels)
labels = ax.get_xticklabels()
labels[0] = ''
ax.set_xticklabels(labels)

plt.savefig('metal_map.pdf')
