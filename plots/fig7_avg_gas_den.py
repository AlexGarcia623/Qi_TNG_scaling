import matplotlib.pyplot as plt
import glob
import os
import sys
import time
import h5py
from scipy.interpolate import UnivariateSpline
from matplotlib.colors import LogNorm
import numpy as np
###################################################################
data_dir = '../avgden/'

h = 0.6774
nbins = 10
log_min = -8.5
log_max = -7.0
do_central = False

file = h5py.File(data_dir + 'spx_data.hdf5')

spx_mass = np.array(file['tot/spx_stellar_mass'])
spx_sfr = np.array(file['tot/spx_sfr'])
spx_gm = np.array(file['tot/spx_gas_mass'])
spx_rad = np.array(file['tot/spx_radius'])
spx_ags = np.array(file['tot/spx_ags'])
spx_agd = np.array(file['tot/spx_agd'])
spx_fg = np.array(file['tot/spx_fg'])

sel_gm, sel_agd, sel_ags, sel_fg, sel_sm = [], [], [], [], []

for i in xrange(len(spx_mass)):
        if (np.log10(spx_agd[i]) < log_max) and (np.log10(spx_agd[i]) > log_min):
                sel_gm.append(spx_gm[i])
                sel_ags.append(spx_ags[i])
                sel_agd.append(spx_agd[i])
                sel_fg.append(spx_fg[i])
                sel_sm.append(spx_mass[i])


sel_gm = np.log10(np.array(sel_gm)) + 6.0
print len(sel_gm)
sel_agd = np.log10(np.array(sel_agd))
#print len(sel_agd)
#sel_agd = [x for x in sel_agd if x > 6.0]
#print len(sel_agd)
sel_sm = np.log10(np.array(sel_sm)) + 6.0
sel_fg = np.log10(np.array(sel_fg))

fig = plt.figure(figsize=(15.0, 15.0))
ax = fig.add_axes([0.13, 0.10, 0.85, 0.85])

p, x, patches = ax.hist(sel_gm, 50, facecolor='b', alpha=0.7)

x = x[:-1] + (x[1] - x[0]) / 2
f = UnivariateSpline(x, p, s=50)
ax.plot(x, f(x), color='black', lw=5)

'''
plt.hist2d(sel_gm, sel_agd, bins=20, norm=LogNorm(), cmap='coolwarm')

ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\star}) / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2}$', fontsize=40)
ax.set_ylabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{gas}}) / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2}$', fontsize=40)

ax.tick_params(labelsize=30)
#ax.set_xlim([6.0, 9.0])
#ax.set_ylim([-4.0, 0.0])

plt.savefig('smgm_test.pdf')
#exit()
'''
##############################################
gal_bin = [9.0, 9.5, 10.0, 10.5, 11.0, 15.0]

for imass in range(len(gal_bin) - 1):
        file = h5py.File(data_dir + 'test_data2.hdf5', 'r')
        name = 'mass' + str(int(gal_bin[imass] * 10)).zfill(3)
        mbin_data = file[name]

        spx_mass = np.array(mbin_data['spx_sm'])
        spx_sfr = np.array(mbin_data['spx_sfr'])
        spx_gm = np.array(mbin_data['spx_gm'])
        spx_rad = np.array(mbin_data['spx_rad'])
        spx_agd = np.array(mbin_data['spx_agd'])
        spx_ags = np.array(mbin_data['spx_ags'])
        spx_fg = np.array(mbin_data['spx_fg'])

        file.close()

        sel_gm, sel_agd, sel_ags, sel_fg, sel_sm = [], [], [], [], []
        for i in xrange(len(spx_mass)):
                if (np.log10(spx_agd[i]) < log_max) and (np.log10(spx_agd[i]) > log_min):
                        sel_gm.append(spx_gm[i])
                        sel_ags.append(spx_ags[i])
                        sel_agd.append(spx_agd[i])
                        sel_fg.append(spx_fg[i])
                        sel_sm.append(spx_mass[i])

        sel_gm = np.log10(np.array(sel_gm)) + 6.0
        sel_agd = np.log10(np.array(sel_agd))
        print len(sel_agd)
        #sel_agd = [x for x in sel_agd if x > 6.0]
        #print len(sel_agd)

        p, x, patches = ax.hist(sel_gm, 30, facecolor='b', alpha=0.7)
        x = x[:-1] + (x[1] - x[0]) / 2
        f = UnivariateSpline(x, p, s=30)
        ax.plot(x, f(x), label= str(gal_bin[imass]) + ' - ' + str(gal_bin[imass+1]), lw=3)

plt.legend(fancybox=True, framealpha=0.5, fontsize=30)
ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{gas}})\,/\,\mathrm{M}_{\odot}\,\mathrm{kpc}^{-2}$', fontsize=40)
ax.set_ylabel('Count', fontsize=40)
ax.tick_params(labelsize=40)

plt.savefig('agd.pdf')

