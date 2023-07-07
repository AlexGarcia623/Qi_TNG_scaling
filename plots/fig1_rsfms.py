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

save_dir = '../'
fig_name = save_dir + 'hist2d_sfr.pdf'

###################################################################
#data_dir = '../fixed/'
#data_dir = '/orange/paul.torrey/jqi/Project/TNG_analysis/grid_maps/L35n2160TNG/'
# data_dir = '../L35n2160TNG/'

data_dir  = '/orange/paul.torrey/jqi/TNG_resolved_scaling_relation/L35n2160TNG/'
data_dir1 = '/orange/paul.torrey/jqi/TNG_resolved_scaling_relation/L75n1820FP/'

# Load in Spaxel Data
file = h5py.File(data_dir + 'spx_data.hdf5', 'r')

grp = file['cen']

spx_mass  = np.array(grp['spx_stellar_mass'])
spx_sfr   = np.array(grp['spx_sfr'])
spx_metal = np.array(grp['spx_metallicity'])
spx_gm    = np.array(grp['spx_gas_mass'])

spx_mass  = np.log10(spx_mass) + 6.0
spx_sfr   = np.log10(spx_sfr)
spx_metal = np.log10(spx_metal * 0.0127 * 0.35 / 16) + 12.0
spx_gm    = np.log10(spx_gm) + 6.0

print 'lowest 3:' , np.sort(spx_metal)[:10]

file.close()

# Load in Median Data for Data Set
mass, mass1, mass2, sfr, sfr1, sfr2, z, z1, z2 = \
        np.loadtxt(data_dir + 'res/res_cen10.txt', unpack=True, usecols=(0,1,2,3,4,5,6,7,8))

mass  = np.log10(mass) + 6.0
mass1 = np.log10(mass1) + 6.0
mass2 = np.log10(mass2) + 6.0
sfr   = np.log10(sfr)
sfr1  = np.log10(sfr1)
sfr2  = np.log10(sfr2)
z     = np.log10(z)
z1    = np.log10(z1)
z2    = np.log10(z2)

# Load in Comparison Data
# TF18
tf18_mass,   tf18_sfr   = np.loadtxt('../tf18.txt', usecols=(0, 1), unpack=True)
tf18_mass_l, tf18_sfr_l = np.loadtxt('../tf181.txt', usecols=(0, 1), unpack=True)
tf18_mass_h, tf18_sfr_h = np.loadtxt('../tf182.txt', usecols=(0, 1), unpack=True)

tf18_sfr   += 6.0
tf18_sfr_l += 6.0
tf18_sfr_h += 6.0
tf18_mass  += 6.0

# MaNGA
manga_mass, manga_sfr = np.loadtxt('../manga.txt', usecols=(0, 1), unpack=True)
manga_sfr  += 6.0
manga_mass += 6.0

#manga_mass = [7.0, 9.0]
#manga_sfr = [7.0 * 0.715 - , 9.0 * 0.715 - ]

# CALIFA
#califa_mass, califa_sfr = np.loadtxt('../califa.txt', usecols=(0, 1), unpack=True)
califa_mass = [7.0, 9.0]
califa_sfr  = [7.0 * 0.72 - 7.95, 9.0 * 0.72 - 7.95]
# SAMI

# Illustris
# data_dir1 = '../L75n1820FP/'
il_mass, il_mass1, il_mass2, il_sfr, il_sfr1, il_sfr2, il_z, il_z1, il_z2 = \
	np.loadtxt(data_dir1 + 'res/res10.txt', unpack=True, usecols=(0,1,2,3,4,5,6,7,8))

il_mass  = np.log10(il_mass) + 6.0
il_mass1 = np.log10(il_mass1) + 6.0
il_mass2 = np.log10(il_mass2) + 6.0
il_sfr   = np.log10(il_sfr)
il_sfr1  = np.log10(il_sfr1)
il_sfr2  = np.log10(il_sfr2)

# Plot Figure
fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.10, 0.10, 0.85, 0.85])

plt.hist2d(spx_mass, spx_sfr, bins=20, norm=LogNorm(), cmap='coolwarm')

ax.plot(mass, sfr, 'black', lw=5, label='IllustrisTNG')

if False:
        ax.plot(il_mass, il_sfr, 'purple', lw=5, label='Illustris')

if True:
	ax.plot(tf18_mass,   tf18_sfr,   color='red', lw=5, label='EAGLE')
	#ax.plot(tf18_mass_l, tf18_sfr_l, color='red', ls='--', lw=3)
	#ax.plot(tf18_mass_h, tf18_sfr_h, color='red', ls='--', lw=3)
if True:
	ax.plot(manga_mass, manga_sfr, color='grey', lw=5, label='MaNGA')
	ax.plot(califa_mass, califa_sfr, color='blue', lw=5, label='CALIFA')

ax.plot(mass, sfr, 'black', lw=5)
ax.plot(mass, sfr1, 'black', ls='--', lw=3)
ax.plot(mass, sfr2, 'black', ls='--', lw=3)

ax.set_xlabel('$\mathrm{log}(\Sigma_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=40)
ax.set_ylabel('$\mathrm{log}(\Sigma_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1})$', fontsize=40)
ax.tick_params(labelsize=40)
cbar = plt.colorbar(pad=0.03)
cbar.set_label('Count', fontsize=50, rotation=90, labelpad=1.1)
cbar.ax.tick_params(labelsize=40)
ax.legend(fancybox=True, framealpha=0.5, fontsize=30)
ax.set_xlim([7.0, 11.2])
ax.set_ylim([-3.9, 1.2])
plt.savefig(fig_name)

plt.clf()

x = spx_mass
y = spx_sfr
z = spx_metal
bins = 40

print z

sums  , xbins, ybins = np.histogram2d(x, y, bins=bins, weights=z)
counts,     _,     _ = np.histogram2d(x, y, bins=(xbins,ybins))

sums   = np.transpose(sums)
counts = np.transpose(counts)

with np.errstate(divide='ignore', invalid='ignore'):  # suppress possible divide-by-zero warnings
    img = plt.pcolor(xbins, ybins, sums / counts, vmin=8.5,vmax=8.8)

ax = plt.gca()

ax.set_xlabel('$\mathrm{log}(\Sigma_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=40)
ax.set_ylabel('$\mathrm{log}(\Sigma_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1})$', fontsize=40)
ax.tick_params(labelsize=40)
cbar = plt.colorbar(pad=0.03)
cbar.set_label(r'$\log{\rm O/H} + 12$', fontsize=50, rotation=90, labelpad=1.1)
cbar.ax.tick_params(labelsize=40)
ax.set_xlim([7.0, 11.2])
ax.set_ylim([-3.9, 1.2])

plt.tight_layout()

plt.savefig(save_dir + 'rFMR.pdf')