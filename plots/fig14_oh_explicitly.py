import matplotlib
matplotlib.use('agg')

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
import h5py
from matplotlib.colors import LogNorm

h = 0.6774
fig_name = 'hist2d_oh.pdf'

###################################################################
data_dir = '../o_over_h/'

# Load in Spaxel Data
file = h5py.File(data_dir + 'spx_data.hdf5', 'r')

grp = file['tot']

spx_mass = np.array(grp['spx_stellar_mass'])
spx_metal = np.array(grp['spx_metallicity'])

spx_mass = np.log10(spx_mass) + 6.0
spx_metal = np.log10(spx_metal) + 12.0
file.close()

# Load in Median Data for Data Set
mass, mass1, mass2, z, z1, z2 = \
        np.loadtxt(data_dir + 'res/res10.txt', unpack=True, usecols=(0,1,2,3,4,5))

mass = np.log10(mass) + 6.0
mass1 = np.log10(mass1) + 6.0
mass2 = np.log10(mass2) + 6.0
z = np.log10(z) + 12.0
z1 = np.log10(z1) + 12.0
z2 = np.log10(z2) + 12.0

###############################################################
tf18_mass, tf18_z = np.loadtxt('tf18_z.txt', usecols=(0, 1), unpack=True)
califa_mass, califa_z = np.loadtxt('califa_z.txt', usecols=(0, 1), unpack=True)
manga_mass, manga_z = np.loadtxt('manga_z.txt', usecols=(0, 1), unpack=True)
il_mass, il_z = np.loadtxt('../L75n1820FP/res/res10.txt', usecols=(0, 6), unpack=True)
il_mass = np.log10(il_mass) + 6.0
il_z = np.log10(il_z * 0.0127 * 0.35 / 16 / 0.75) + 12.0
###############################################################
###############################################################
fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.15, 0.10, 0.85, 0.85])

plt.hist2d(spx_mass, spx_metal, range=[[7.0, 9.5], [8.0, 9.25]], bins=20, norm=LogNorm(), cmap='coolwarm')

ax.plot(mass, z, 'black', lw=5, label='IllustrisTNG')
ax.plot(il_mass, il_z, color='purple', lw=5, label='Illustris')

ax.plot(tf18_mass, tf18_z, color='red', lw=5, label='EAGLE')
ax.plot(manga_mass, manga_z, color='grey', lw=5, label='MaNGA')
ax.plot(califa_mass, califa_z, color='blue', lw=5, label='CALIFA')

ax.plot(mass, z, 'black', lw=5)
ax.plot(mass, z1, 'black', ls='--', lw=3)
ax.plot(mass, z2, 'black', ls='--', lw=3)

ax.set_xlabel('$\mathrm{log}(\Sigma_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=40)
ax.set_ylabel('$12 + \mathrm{log}(\mathrm{O}/\mathrm{H})$', fontsize=40)
ax.set_ylim([8.00, 9.60])
ax.set_xlim([7.01, 9.50])
ax.tick_params(labelsize=40)
cbar = plt.colorbar(pad=0.03)
cbar.set_label('Count', fontsize=50, rotation=90, labelpad=1.1)
cbar.ax.tick_params(labelsize=40)
ax.legend(fancybox=True, framealpha=0.5, fontsize=30, loc='upper left')

plt.savefig(fig_name)
