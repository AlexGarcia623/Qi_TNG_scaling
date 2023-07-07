import matplotlib
matplotlib.use('agg')

import sys
utils   = '/home/alexgarcia/torreylabtools/Python/'

sys.path.insert(2,utils)

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
from util import cosmo_tools
import h5py
from matplotlib.colors import LogNorm
####################################################################################
# data_dir = '/orange/paul.torrey/jqi/fixed/'
data_dir  = '/orange/paul.torrey/jqi/TNG_resolved_scaling_relation/fixed_x/'

# Load in Spaxel Data
file = h5py.File(data_dir + 'spx_data.hdf5', 'r')

count = 0
for j in file:
    count += 1

print count

grp = file['tot']

spx_sm = np.array(grp['spx_stellar_mass'])
spx_sfr = np.array(grp['spx_sfr'])
spx_metal = np.array(grp['spx_metallicity'])
spx_gm = np.array(grp['spx_gas_mass'])

spx_fg = np.log10( spx_gm / (spx_gm + spx_sm) )

spx_sm = np.log10(spx_sm) + 6.0
spx_sfr = np.log10(spx_sfr)
spx_gm = np.log10(spx_gm) + 6.0

file.close()

print len(spx_fg)

###################################################################################

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.85, 0.85])

plt.hist2d(spx_sm, spx_fg, bins=20, norm=LogNorm(), cmap='coolwarm')

# Load in Median Data for Data Set
mass, mass1, mass2, sfr, sfr1, sfr2 = \
        np.loadtxt(data_dir + 'res/res_ks210.txt', unpack=True, usecols=(0,1,2,3,4,5))

mass = np.log10(mass) + 6.0
mass1 = np.log10(mass1) + 6.0
mass2 = np.log10(mass2) + 6.0
sfr = np.log10(sfr)
sfr1 = np.log10(sfr1)
sfr2 = np.log10(sfr2)

print mass, sfr

ax.plot(mass, sfr, 'black', lw=5)
#ax.plot(mass, sfr1, 'black', ls='--', lw=3)
#ax.plot(mass, sfr2, 'black', ls='--', lw=3)

file_list = sorted(glob.glob(data_dir + 'res/ks2_bin???.txt'))
print file_list

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0, 15.0]
colors = ['red', 'orange', 'green', 'blue']
for i in xrange(len(file_list)):
        m, m1, m2, s, s1, s2 = \
                np.loadtxt(file_list[i], unpack=True, usecols=(0,1,2,3,4,5))
        m = np.log10(m) + 6.0
        s = np.log10(s)
        if (i < len(file_list) - 1):
                label = '$ 10^{' + str(mass_bin[i]) + '}\,<\,M_{\star}\,<\,10^{' + str(mass_bin[i+1]) + '}\,\mathrm{M}_{\odot}$'
                ax.plot(m, s, color=colors[i], label= label, lw=5) #str(mass_bin[i]) + ' - ' + str(mass_bin[i+1]), lw=5)


ax.set_xlabel('$\mathrm{log}(\Sigma_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=40)
ax.set_ylabel('$\mathrm{log}(f_{\mathrm{gas}})$', fontsize=40)
ax.tick_params(labelsize=40)
cbar = plt.colorbar(pad=0.03)
cbar.set_label('Count', fontsize=50, rotation=90, labelpad=1.1)
cbar.ax.tick_params(labelsize=40)
#ax.legend(fancybox=True, framealpha=0.5, fontsize=30)
ax.set_xlim([7.1, 11.2])
ax.set_ylim([-4.1, -0.1])
plt.savefig('test.png')

