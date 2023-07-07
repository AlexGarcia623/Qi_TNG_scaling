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
data_dir = '../L35n2160TNG/'
h = 0.6774

mass, mass1, mass2, sfr, sfr1, sfr2, z, z1, z2 = \
        np.loadtxt(data_dir + 'res/res_cen10.txt', unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
mass = np.log10(mass) + 6.0
sfr = np.log10(sfr)
z = np.log10(z)

file_list = sorted(glob.glob(data_dir + 'res/test_res_bin???_cen.txt'))
print file_list

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0, 15.0]
colors = ['red', 'orange', 'green', 'blue']

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.85, 0.85])
ax.plot(mass, sfr, 'black', lw=5)

for i in xrange(len(file_list)):
        m, m1, m2, s, s1, s2, z, z1, z2 = \
                np.loadtxt(file_list[i], unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
        m = np.log10(m) + 6.0
        s = np.log10(s)
        z = np.log10(z)
        if (i < len(file_list) - 1):
		label = '$ 10^{' + str(mass_bin[i]) + '}\,<\,M_{\star}\,<\,10^{' + str(mass_bin[i+1]) + '}\,\mathrm{M}_{\odot}$'
                ax.plot(m, s, color=colors[i], label= label, lw=5) #str(mass_bin[i]) + ' - ' + str(mass_bin[i+1]), lw=5)

ax.set_xlabel('$\mathrm{log}(\Sigma_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=40)
ax.set_ylabel('$\mathrm{log}(\Sigma_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1}) $', fontsize=40)
ax.set_xlim([7.01, 9.7])
ax.set_ylim([-4.1, -0.4])
ax.tick_params(labelsize=40)
plt.legend(fontsize=40)
plt.savefig('bin_rsfms.pdf')
