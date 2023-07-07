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

###################################################################
data_dir = '../fixed/'
h = 0.6774

mass, mass1, mass2, sfr, sfr1, sfr2, z, z1, z2 = \
        np.loadtxt(data_dir + 'res/res_cen10.txt', unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
mass = np.log10(mass) + 6.0
sfr = np.log10(sfr)
z = np.log10(z * 0.0127 * 0.35 / 16) + 12.0
mass1 = np.log10(mass1) + 6.0
mass2 = np.log10(mass2) + 6.0
z1 = np.log10(z1 * 0.0127 * 0.35 / 16) + 12.0
z2 = np.log10(z2 * 0.0127 * 0.35 / 16) + 12.0

file_list = sorted(glob.glob(data_dir + 'res/test_res_bin???_cen.txt'))
print file_list

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0, 15.0]
colors = ['red', 'orange', 'green', 'blue']

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.85, 0.85])
ax.plot(mass, z, 'black', lw=5)
ax.plot(mass1, z1, 'black', ls='--', lw=3)
ax.plot(mass2, z2, 'black', ls='--', lw=3)

for i in xrange(len(file_list)):
        m, m1, m2, s, s1, s2, z, z1, z2 = \
                np.loadtxt(file_list[i], unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
        m = np.log10(m) + 6.0
        s = np.log10(s * h)
        z = np.log10(z * 0.0127 * 0.35 / 16) + 12.0
        if (i < len(file_list) - 1):
                label = '$ 10^{' + str(mass_bin[i]) + '}\,<\,M_{\star}\,<\,10^{' + str(mass_bin[i+1]) + '}\,\mathrm{M}_{\odot}$'
                ax.plot(m, z, color=colors[i], label=label, lw=5)

ax.set_xlabel('$\mathrm{log}(\Sigma_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=40)
ax.set_ylabel('$12 + \mathrm{log}(\mathrm{O}/\mathrm{H})$', fontsize=40)
ax.set_ylim([8.00, 9.60])
ax.set_xlim([7.01, 9.50])
ax.tick_params(labelsize=40)
plt.legend(fontsize=40)
plt.savefig('bin_rmzr.pdf')
