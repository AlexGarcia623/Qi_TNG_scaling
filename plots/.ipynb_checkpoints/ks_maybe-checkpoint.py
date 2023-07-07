import matplotlib as mpl
mpl.use('agg')

mpl.rcParams['font.size']=50

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import time
import h5py
from matplotlib.colors import LogNorm
h = 0.6774

###################################################################
# data_dir = '../fixed/'

data_dir = '/orange/paul.torrey/jqi/TNG_resolved_scaling_relation/fixed_x/'
save_dir = '../'

STARS_OR_GAS = 'gas'.upper()

file = h5py.File(data_dir + 'test_data2.hdf5', 'r')

file_list = sorted(glob.glob(data_dir + 'res/test_res_bin???.txt'))
names = ['mass090', 'mass095', 'mass100', 'mass105', 'mass110']
print file_list
print names

titles = ['$10^{9.0}\,<\,\mathrm{M}_{\star}\,<\,10^{9.5}\,\mathrm{M}_{\odot}$', \
        '$10^{9.5}\,<\,\mathrm{M}_{\star}\,<\,10^{10.0}\,\mathrm{M}_{\odot}$', \
        '$10^{10.0}\,<\,\mathrm{M}_{\star}\,<\,10^{10.5}\,\mathrm{M}_{\odot}$', \
        '$10^{10.5}\,<\,\mathrm{M}_{\star}\,<\,10^{11.0}\,\mathrm{M}_{\odot}$', \
        '$10^{11.0}\,<\,\mathrm{M}_{\star}\,<\,10^{15.0}\,\mathrm{M}_{\odot}$']

fig, axs = plt.subplots(2,2, figsize=(30.0,30.0), sharex=True, sharey=True)
axs = axs.flatten()


fig_name = save_dir + 'ks_maybe.pdf'

for k in xrange(len(names) - 1):
    grp = file[names[k]]
    
    spx_sm    = np.log10( np.array(grp['spx_stellar_mass'])) + 6.0
    spx_sfr   = np.log10( np.array(grp['spx_sfr']) * h     )
    spx_metal = np.log10( np.array(grp['spx_metallicity']) )
    spx_gm    = np.log10( np.array(grp['spx_gas_mass'])    ) + 6.0
    spx_rad   =           np.array(grp['spx_radius'])
    
    sfr_l = -4.0
    sfr_h = 0.0
    sm_l  = 7.0
    sm_h  = 10.0

    if (STARS_OR_GAS == 'gas'.upper()):
        sm_l = 6.0
        sm_h = 9.0
    
    nbins  = 20
    sfr_sp = (sfr_h - sfr_l) / nbins
    sm_sp  = (sm_h - sm_l) / nbins
    m2     = np.zeros((nbins, nbins))

    m1 = [[[] for _ in range(nbins) ] for _ in range(nbins)]

    spx_mass = spx_sm
    if (STARS_OR_GAS == 'gas'.upper()):
        spx_mass = spx_gm
    
    for i in xrange(len(spx_mass)):
        if spx_sfr[i] < sfr_l: continue
        y = nbins - 1 - int( (spx_sfr[i] - sfr_l) / sfr_sp )
        x = int( (spx_mass[i] - sm_l) / sm_sp )
        if (x >= 0) and (x < nbins) and (y >= 0) and (y < nbins):
            m1[y][x].append(spx_rad[i])
            m2[y][x] += 1


    m3 = np.zeros((nbins, nbins))
    for i in range(nbins):
        for j in range(nbins):
            if m2[i][j] > 10:
                m3[i][j] = np.median(m1[i][j])
            else:
                m3[i][j] = None

    mass, mass1, mass2, sfr, sfr1, sfr2, z, z1, z2 = \
            np.loadtxt(file_list[k], unpack=True, usecols=(0,1,2,3,4,5,6,7,8))
    mass = np.log10(mass) + 6.0
    sfr  = np.log10(sfr * h)
    sfr1 = np.log10(sfr1 * h)
    sfr2 = np.log10(sfr2 * h)
    
    # fig = plt.figure(figsize=(15.0, 15.0))
    # ax = fig.add_axes([0.16, 0.10, 0.82, 0.82])

    ax = axs[k]

    mpb = ax.imshow(m3, extent=(sm_l, sm_h, sfr_l, sfr_h), interpolation='none', vmin=0, vmax=12, cmap='BrBG', aspect='auto')
    # ax.plot(mass, sfr, 'black', lw=5, label='IllustrisTNG')
    # ax.plot(mass, sfr1, 'black', ls='--', lw=3)
    # ax.plot(mass, sfr2, 'black', ls='--', lw=3)
    
    if (k == 0 or k==2):
        ax.set_ylabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{SFR}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2} \, \mathrm{yr}^{-1})$', fontsize=50)
    if (k == 2 or k==3):
        ax.set_xlabel('$\mathrm{log}(\mathrm{\Sigma}_{\mathrm{\star}} / \mathrm{M}_{\odot}\, \mathrm{kpc}^{-2})$', fontsize=50)

    ax.tick_params(labelsize=50)
    # ax.set_xlim([7.01, 9.99])
    ax.set_ylim([-4.0, -0.01])

    ax.text( 0.45,0.9, titles[k] ,transform=ax.transAxes, fontsize=50 )

ax_cbar = fig.add_axes([0.2, 1., 0.6, 0.02])

cb = plt.colorbar(mpb, cax=ax_cbar,orientation='horizontal',aspect=10,shrink=0.5)
# cb.set_label(r'${\rm Radius~(kpc)}$',fontsize=50,labelpad=100)

cb.ax.text( 0.5,1.25, r'${\rm Radius~(kpc)}$' , transform=cb.ax.transAxes, fontsize=50, ha='center' )

# cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')

cb.ax.tick_params(labelsize=50)

plt.tight_layout()

plt.subplots_adjust(wspace=0.01, hspace=0.01)
plt.savefig(fig_name,bbox_inches='tight')

file.close()