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
data_dir = '../appendix/only_sfr_gas/'

#file = h5py.File( data_dir + 'spx_data.hdf5', 'r')

file_list = sorted(glob.glob(data_dir + 'res/test_res_bin???.txt'))
print file_list

mass_bin = [9.0, 9.5, 10.0, 10.5, 11.0, 15.0]
colors = ['red', 'orange', 'green', 'blue']

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.85, 0.85])

for i in xrange(len(file_list) - 1):
        sm, gm = \
                np.loadtxt(file_list[i], unpack=True, usecols=(0,9))
        sm = np.log10(sm) + 6.0
        gm = np.log10(gm) + 6.0

        if (i < len(file_list) - 1):
                label = '$ 10^{' + str(mass_bin[i]) + '}\,<\,M_{\star}\,<\,10^{' + str(mass_bin[i+1]) + '}\,\mathrm{M}_{\odot}$'		
                ax.plot(sm, gm, color=colors[i], label=label, lw=5)

ax.set_xlabel('$\mathrm{log}(\Sigma_{\star} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=40)
ax.set_ylabel('$\mathrm{log}(\Sigma_{\mathrm{gas,\,sf}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1})$', fontsize=40)
ax.tick_params(labelsize=40)
plt.legend(fontsize=40)
ax.set_xlim([7.01, 9.7])
ax.set_ylim([6.01, 9.7])
plt.savefig('only_sfr_gas.pdf')
